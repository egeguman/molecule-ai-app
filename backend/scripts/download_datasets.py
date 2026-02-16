"""Download all 19 ADMET benchmark datasets from Therapeutics Data Commons (TDC).

Usage:
    cd backend
    python -m scripts.download_datasets
"""

import os
import sys

# Add backend root to path for config imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from config.settings import ADMET_DATASETS_DIR, ADMET_ENDPOINTS


def download_all():
    """Download all ADMET datasets using PyTDC."""
    from tdc.single_pred import ADME, Tox

    tdc_classes = {"ADME": ADME, "Tox": Tox}

    os.makedirs(ADMET_DATASETS_DIR, exist_ok=True)

    for name, meta in ADMET_ENDPOINTS.items():
        tdc_class = tdc_classes[meta["tdc_class"]]
        tdc_name = meta["tdc_name"]

        print(f"\nDownloading {name} ({meta['tdc_class']}: {tdc_name})...")

        try:
            data = tdc_class(name=tdc_name)
            split = data.get_split()

            dataset_dir = ADMET_DATASETS_DIR / name
            os.makedirs(dataset_dir, exist_ok=True)

            for subset_name, subset_df in split.items():
                output_path = dataset_dir / f"{subset_name}.csv"
                subset_df.to_csv(output_path, index=False)
                print(f"  {subset_name}: {len(subset_df)} samples -> {output_path}")

        except Exception as e:
            print(f"  ERROR downloading {name}: {e}")
            continue

    print(f"\nDone. Datasets saved to {ADMET_DATASETS_DIR}")


if __name__ == "__main__":
    download_all()
