# Claude Code Prompt --- Integrate Flask + Chemprop ADMET + OpenAI GPT-5.2 + Word Report

You are Claude Code running inside my repository. Your task is to
integrate a Flask backend into my existing functional frontend. The
backend must:

1.  Accept a SMILES string
2.  Download and prepare ADMET benchmark datasets
3.  Fine-tune Chemprop models on each ADMET endpoint
4.  Run Chemprop ADMET predictions using the fine-tuned models
5.  Store results in CSV
6.  Send structured results to OpenAI GPT-5.2
7.  Generate and return a Word (.docx) report

------------------------------------------------------------------------

## Backend Requirements

### Endpoints

-   POST /api/predict\
    Input: { "smiles": "CCO" }

-   POST /api/train\
    Triggers fine-tuning pipeline for all ADMET endpoints.\
    Input (optional): { "endpoints": ["ames", "caco2"], "epochs": 50 }

-   GET /api/models/status\
    Returns training status and test metrics for all endpoints.

-   GET /api/health\
    Returns: { "status": "ok" }

------------------------------------------------------------------------

## ADMET Dataset Sourcing

Before training, download and prepare benchmark ADMET datasets. Use the
Therapeutics Data Commons (TDC) Python package as the primary source:

    pip install PyTDC

### Required Datasets (grouped by ADMET category)

**Absorption:**
-   Caco-2 permeability (TDC: `Caco2_Wang`)
-   HIA — Human Intestinal Absorption (TDC: `HIA_Hou`)
-   Pgp Inhibition (TDC: `Pgp_Broccatelli`)
-   Bioavailability (TDC: `Bioavailability_Ma`)

**Distribution:**
-   BBB — Blood-Brain Barrier penetration (TDC: `BBB_Martins`)
-   PPB — Plasma Protein Binding (TDC: `PPBR_AZ`)
-   VDss — Volume of Distribution (TDC: `VDss_Lombardo`)

**Metabolism:**
-   CYP2D6 Inhibition (TDC: `CYP2D6_Veith`)
-   CYP3A4 Inhibition (TDC: `CYP3A4_Veith`)
-   CYP2C9 Inhibition (TDC: `CYP2C9_Veith`)
-   CYP2D6 Substrate (TDC: `CYP2D6_Substrate_CarbonMangels`)
-   CYP3A4 Substrate (TDC: `CYP3A4_Substrate_CarbonMangels`)

**Excretion:**
-   Half-life (TDC: `Half_Life_Obach`)
-   Clearance — Hepatocyte (TDC: `Clearance_Hepatocyte_AZ`)
-   Clearance — Microsome (TDC: `Clearance_Microsome_AZ`)

**Toxicity:**
-   hERG Cardiotoxicity (TDC: `hERG`)
-   AMES Mutagenicity (TDC: `AMES`)
-   DILI — Drug-Induced Liver Injury (TDC: `DILI`)
-   LD50 Acute Toxicity (TDC: `LD50_Zhu`)

### Download Script

Create `backend/scripts/download_datasets.py`:

```python
from tdc.benchmark_group import admet_group
group = admet_group(path='backend/data/admet_datasets/')
# Iterate and export each benchmark to CSV
for task_name in group.dataset_names:
    benchmark = group.get(task_name)
    train, valid = benchmark['train_val'], benchmark['test']
    train.to_csv(f'backend/data/admet_datasets/{task_name}_train.csv', index=False)
    valid.to_csv(f'backend/data/admet_datasets/{task_name}_test.csv', index=False)
```

Store all raw datasets under `backend/data/admet_datasets/`. Each CSV must
have at minimum a `Drug` (SMILES) column and a `Y` (label/value) column.

------------------------------------------------------------------------

## Chemprop Fine-Tuning

Use Chemprop v2 (`chemprop>=2.0.0`) built on PyTorch Lightning.

    pip install chemprop

### Training Pipeline

Create `backend/scripts/train_models.py` that trains one model per ADMET
endpoint:

1.  **For each dataset CSV** in `backend/data/admet_datasets/`:
    -   Detect task type: classification (binary Y) vs regression (continuous Y)
    -   Split: use the TDC-provided scaffold split (already separated into
        train/test); further split train into train/val (80/20 scaffold split)

2.  **Training command per endpoint** (CLI equivalent):

    ```bash
    # Classification example (e.g., AMES)
    chemprop train \
        --data-path backend/data/admet_datasets/AMES_train.csv \
        --smiles-column Drug --target-columns Y \
        --task-type classification \
        --output-dir backend/models/ames/ \
        --epochs 50 --batch-size 64 --ensemble-size 5

    # Regression example (e.g., Caco2)
    chemprop train \
        --data-path backend/data/admet_datasets/Caco2_Wang_train.csv \
        --smiles-column Drug --target-columns Y \
        --task-type regression \
        --output-dir backend/models/caco2/ \
        --epochs 50 --batch-size 64 --ensemble-size 5
    ```

3.  **Ensemble**: Train 5-model ensembles per endpoint (`--ensemble-size 5`)
    to get uncertainty estimates. Average predictions across ensemble members.

4.  **Extra features** (optional but recommended): Append RDKit 2D descriptors
    alongside the learned D-MPNN features:

    ```bash
    chemprop train ... --molecule-featurizers rdkit_2d_normalized
    ```

5.  **Hyperparameter optimization** (optional): Use Chemprop's built-in
    Optuna-based hyperparameter search:

    ```bash
    chemprop hpopt \
        --data-path backend/data/admet_datasets/AMES_train.csv \
        --smiles-column Drug --target-columns Y \
        --task-type classification \
        --hpopt-save-dir backend/models/ames_hpopt/ \
        --num-iters 30
    ```

### Model Storage

Save all trained checkpoints to `backend/models/<endpoint_name>/`. Each
folder should contain:
-   `model_0/` through `model_4/` (ensemble members)
-   `training_args.json` (reproducibility)
-   `test_metrics.json` (performance on held-out TDC test set)

### Evaluation

After training, evaluate every model on the TDC test split and log:
-   AUROC for classification tasks
-   RMSE and MAE for regression tasks
-   Save to `backend/models/<endpoint_name>/test_metrics.json`

------------------------------------------------------------------------

## Chemprop Inference Integration

-   Load all trained Chemprop checkpoints from `backend/models/`
-   Validate SMILES using RDKit (`Chem.MolFromSmiles()` — reject if None)
-   Run inference across all ADMET endpoints in a single request
-   For ensemble models, return mean prediction and standard deviation
-   Return structured ADMET dictionary:

```json
{
  "absorption": {
    "caco2_permeability": {"value": -5.12, "unit": "log cm/s", "std": 0.08},
    "hia": {"value": 0.94, "unit": "probability", "std": 0.02},
    "pgp_inhibitor": {"value": 0.15, "unit": "probability", "std": 0.04},
    "bioavailability": {"value": 0.82, "unit": "probability", "std": 0.03}
  },
  "distribution": {
    "bbb_penetration": {"value": 0.88, "unit": "probability", "std": 0.05},
    "ppb": {"value": 85.2, "unit": "percent", "std": 2.1},
    "vdss": {"value": 0.45, "unit": "log L/kg", "std": 0.12}
  },
  "metabolism": {
    "cyp2d6_inhibition": {"value": 0.08, "unit": "probability", "std": 0.03},
    "cyp3a4_inhibition": {"value": 0.22, "unit": "probability", "std": 0.06},
    "cyp2c9_inhibition": {"value": 0.11, "unit": "probability", "std": 0.04},
    "cyp2d6_substrate": {"value": 0.35, "unit": "probability", "std": 0.07},
    "cyp3a4_substrate": {"value": 0.67, "unit": "probability", "std": 0.05}
  },
  "excretion": {
    "half_life": {"value": 3.2, "unit": "hours", "std": 0.5},
    "clearance_hepatocyte": {"value": 42.1, "unit": "uL/min/10^6 cells", "std": 3.8},
    "clearance_microsome": {"value": 55.3, "unit": "uL/min/mg", "std": 4.2}
  },
  "toxicity": {
    "herg_inhibition": {"value": 0.05, "unit": "probability", "std": 0.02},
    "ames_mutagenicity": {"value": 0.12, "unit": "probability", "std": 0.03},
    "dili": {"value": 0.08, "unit": "probability", "std": 0.02},
    "ld50": {"value": 2.85, "unit": "log mol/kg", "std": 0.15}
  }
}
```

------------------------------------------------------------------------

## CSV Logging

Append to backend/data/predictions.csv with:

-   timestamp_iso
-   request_id
-   smiles
-   predictions_json
-   ensemble_std_json
-   gpt_summary_json
-   latency metrics
-   status
-   error_message

------------------------------------------------------------------------

## OpenAI GPT-5.2 Integration

Use environment variables:

-   OPENAI_API_KEY
-   OPENAI_MODEL (default: gpt-5.2)

Use structured JSON output schema:

```json
{
  "executive_summary": "string",
  "high_risk_flags": ["string"],
  "admet_interpretation": [
    {
      "property": "string",
      "value": "number|string",
      "std": "number",
      "interpretation": "string",
      "confidence_notes": "string"
    }
  ],
  "recommended_next_steps": ["string"],
  "disclaimer": "string"
}
```

------------------------------------------------------------------------

## Word Report Generation

Use python-docx.

Include:

1.  Title
2.  Metadata (request_id, timestamp, SMILES, 2D molecule image via RDKit)
3.  ADMET results table (with values, units, and uncertainty)
4.  Model confidence section (ensemble std per endpoint)
5.  GPT interpretation sections
6.  Disclaimer

Return as downloadable:

Content-Type:
application/vnd.openxmlformats-officedocument.wordprocessingml.document
Content-Disposition: attachment

------------------------------------------------------------------------

## Project Structure

```
backend/
    app.py
    services/
        chemprop_service.py     # inference + model loading
        training_service.py     # fine-tuning pipeline
        openai_service.py       # GPT-5.2 integration
        report_service.py       # Word doc generation
    config/
    models/                     # trained checkpoints per endpoint
        ames/
        caco2/
        bbb/
        ...
    data/
        admet_datasets/         # TDC benchmark CSVs
        predictions.csv         # prediction log
    scripts/
        download_datasets.py    # TDC dataset download
        train_models.py         # full training pipeline
    tests/
```

------------------------------------------------------------------------

## Deliverables

-   Flask backend with /predict, /train, /models/status, /health endpoints
-   Dataset download script (TDC ADMET benchmarks)
-   Chemprop fine-tuning pipeline (per-endpoint, ensemble training)
-   Trained model evaluation (AUROC / RMSE on TDC test splits)
-   Chemprop inference service (ensemble predictions with uncertainty)
-   CSV logging (including ensemble std)
-   OpenAI structured summary
-   Word report generation (with molecule image and confidence metrics)
-   Frontend integration
-   Tests
-   README (including instructions for dataset download and model training)
