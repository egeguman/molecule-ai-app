"""Train Chemprop v2 models for all ADMET endpoints.

Trains a 5-model ensemble per endpoint using PyTorch Lightning.

Usage:
    cd backend
    python -m scripts.train_models
    python -m scripts.train_models --endpoints AMES hERG --epochs 50
"""

import argparse
import json
import logging
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import pandas as pd
import pytorch_lightning as pl
from sklearn.metrics import mean_absolute_error, mean_squared_error, roc_auc_score

from config.settings import (
    ADMET_DATASETS_DIR,
    ADMET_ENDPOINTS,
    ENSEMBLE_SIZE,
    MODELS_DIR,
    TRAINING_BATCH_SIZE,
    TRAINING_EPOCHS,
)

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def train_endpoint(endpoint_name: str, meta: dict, epochs: int | None = None):
    """Train a Chemprop ensemble for a single ADMET endpoint."""
    from chemprop.data import MoleculeDataModule, MoleculeDatapoint, MoleculeDataset
    from chemprop.nn import (
        BinaryClassificationFFN,
        BondMessagePassing,
        MoleculeModel,
        RegressionFFN,
    )

    epochs = epochs or TRAINING_EPOCHS

    # Load data
    data_dir = ADMET_DATASETS_DIR / endpoint_name
    train_path = data_dir / "train.csv"
    valid_path = data_dir / "valid.csv"
    test_path = data_dir / "test.csv"

    if not train_path.exists():
        logger.error("Training data not found at %s. Run download_datasets.py first.", train_path)
        return None

    train_df = pd.read_csv(train_path)
    val_df = pd.read_csv(valid_path) if valid_path.exists() else train_df.sample(frac=0.2)
    test_df = pd.read_csv(test_path) if test_path.exists() else None

    # Build Chemprop datapoints
    smiles_col = "Drug"
    target_col = "Y"

    train_data = [
        MoleculeDatapoint(row[smiles_col], [float(row[target_col])])
        for _, row in train_df.iterrows()
        if pd.notna(row[target_col])
    ]
    val_data = [
        MoleculeDatapoint(row[smiles_col], [float(row[target_col])])
        for _, row in val_df.iterrows()
        if pd.notna(row[target_col])
    ]

    output_dir = MODELS_DIR / endpoint_name
    output_dir.mkdir(parents=True, exist_ok=True)

    logger.info(
        "Training %s: %d train, %d val samples, %d epochs, ensemble=%d",
        endpoint_name, len(train_data), len(val_data), epochs, ENSEMBLE_SIZE,
    )

    # Train ensemble
    for i in range(ENSEMBLE_SIZE):
        logger.info("  Ensemble member %d/%d", i + 1, ENSEMBLE_SIZE)

        mp = BondMessagePassing()
        if meta["type"] == "classification":
            ffn = BinaryClassificationFFN()
        else:
            ffn = RegressionFFN()

        model = MoleculeModel(message_passing=mp, output_transform=ffn)

        datamodule = MoleculeDataModule(
            train=MoleculeDataset(train_data),
            val=MoleculeDataset(val_data),
            batch_size=TRAINING_BATCH_SIZE,
            num_workers=0,
        )

        trainer = pl.Trainer(
            max_epochs=epochs,
            accelerator="auto",
            default_root_dir=str(output_dir),
            enable_progress_bar=True,
            callbacks=[
                pl.callbacks.ModelCheckpoint(
                    dirpath=str(output_dir),
                    filename=f"ensemble_{i}",
                    monitor="val_loss",
                    save_top_k=1,
                    mode="min",
                ),
                pl.callbacks.EarlyStopping(
                    monitor="val_loss",
                    patience=5,
                    mode="min",
                ),
            ],
            logger=False,
        )
        trainer.fit(model, datamodule)

    # Evaluate on test set
    if test_df is not None:
        metrics = _evaluate(endpoint_name, meta, test_df, output_dir)
        metrics_path = output_dir / "test_metrics.json"
        with open(metrics_path, "w") as f:
            json.dump(metrics, f, indent=2)
        logger.info("  Test metrics: %s", metrics)
        return metrics

    return None


def _evaluate(endpoint_name: str, meta: dict, test_df: pd.DataFrame, model_dir) -> dict:
    """Evaluate ensemble on test set."""
    import math

    import torch
    from chemprop.data import MoleculeDatapoint, MoleculeDataset, collate_batch
    from chemprop.nn import MoleculeModel

    ckpts = sorted(model_dir.glob("ensemble_*.ckpt"))
    if not ckpts:
        return {"error": "No checkpoints found"}

    test_data = [
        MoleculeDatapoint(row["Drug"], [float(row["Y"])])
        for _, row in test_df.iterrows()
        if pd.notna(row["Y"])
    ]
    dataset = MoleculeDataset(test_data)
    targets = [row["Y"] for _, row in test_df.iterrows() if pd.notna(row["Y"])]

    # Collect ensemble predictions
    all_preds = []
    for ckpt_path in ckpts:
        model = MoleculeModel.load_from_checkpoint(str(ckpt_path))
        model.eval()
        preds = []
        with torch.no_grad():
            for dp in dataset:
                batch = collate_batch([dp])
                pred = model(batch)
                preds.append(pred.item())
        all_preds.append(preds)

    # Average ensemble predictions
    n = len(targets)
    mean_preds = [sum(all_preds[j][i] for j in range(len(ckpts))) / len(ckpts) for i in range(n)]

    metrics = {"endpoint": endpoint_name, "type": meta["type"], "n_test": n}

    if meta["type"] == "classification":
        # Apply sigmoid
        mean_preds_prob = [1 / (1 + math.exp(-p)) for p in mean_preds]
        try:
            metrics["auroc"] = round(roc_auc_score(targets, mean_preds_prob), 4)
        except ValueError:
            metrics["auroc"] = None
    else:
        import numpy as np
        metrics["rmse"] = round(float(np.sqrt(mean_squared_error(targets, mean_preds))), 4)
        metrics["mae"] = round(float(mean_absolute_error(targets, mean_preds)), 4)

    return metrics


def train_all(endpoints: list[str] | None = None, epochs: int | None = None):
    """Train models for all (or selected) endpoints."""
    targets = endpoints or list(ADMET_ENDPOINTS.keys())
    results = {}

    for name in targets:
        if name not in ADMET_ENDPOINTS:
            logger.warning("Unknown endpoint: %s, skipping", name)
            continue

        logger.info("\n%s", "=" * 60)
        logger.info("Training %s (%s)", name, ADMET_ENDPOINTS[name]["type"])
        logger.info("=" * 60)

        try:
            metrics = train_endpoint(name, ADMET_ENDPOINTS[name], epochs=epochs)
            results[name] = {"status": "completed", "metrics": metrics}
        except Exception as e:
            logger.error("Failed to train %s: %s", name, e)
            results[name] = {"status": "failed", "error": str(e)}

    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train Chemprop ADMET models")
    parser.add_argument("--endpoints", nargs="+", help="Specific endpoints to train")
    parser.add_argument("--epochs", type=int, help="Override training epochs")
    args = parser.parse_args()

    train_all(endpoints=args.endpoints, epochs=args.epochs)
