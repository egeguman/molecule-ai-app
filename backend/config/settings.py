import os
from pathlib import Path
from dotenv import load_dotenv

load_dotenv(Path(__file__).resolve().parent.parent / ".env")

BASE_DIR = Path(__file__).resolve().parent.parent
MODELS_DIR = BASE_DIR / "models"
DATA_DIR = BASE_DIR / "data"
ADMET_DATASETS_DIR = DATA_DIR / "admet_datasets"
PREDICTIONS_CSV = DATA_DIR / "predictions.csv"

OPENAI_API_KEY = os.getenv("OPENAI_API_KEY", "")
OPENAI_MODEL = os.getenv("OPENAI_MODEL", "gpt-4o")

# Auth / Database
JWT_SECRET_KEY = os.getenv("JWT_SECRET_KEY", "dev-fallback-secret-change-me")
JWT_ACCESS_TOKEN_EXPIRES_HOURS = int(os.getenv("JWT_ACCESS_TOKEN_EXPIRES_HOURS", "24"))
DATABASE_URL = os.getenv("DATABASE_URL", f"sqlite:///{BASE_DIR / 'moleculeai.db'}")

ENSEMBLE_SIZE = 5
TRAINING_EPOCHS = int(os.getenv("TRAINING_EPOCHS", "30"))
TRAINING_BATCH_SIZE = int(os.getenv("TRAINING_BATCH_SIZE", "64"))

# All 19 ADMET endpoints with metadata
ADMET_ENDPOINTS = {
    # Absorption
    "Caco2_Wang": {
        "type": "regression",
        "unit": "log cm/s",
        "group": "absorption",
        "display_name": "Caco-2 Permeability",
        "tdc_name": "Caco2_Wang",
        "tdc_class": "ADME",
    },
    "HIA_Hou": {
        "type": "classification",
        "unit": "probability",
        "group": "absorption",
        "display_name": "Human Intestinal Absorption",
        "tdc_name": "HIA_Hou",
        "tdc_class": "ADME",
    },
    "Pgp_Broccatelli": {
        "type": "classification",
        "unit": "probability",
        "group": "absorption",
        "display_name": "P-gp Inhibition",
        "tdc_name": "Pgp_Broccatelli",
        "tdc_class": "ADME",
    },
    "Bioavailability_Ma": {
        "type": "classification",
        "unit": "probability",
        "group": "absorption",
        "display_name": "Bioavailability",
        "tdc_name": "Bioavailability_Ma",
        "tdc_class": "ADME",
    },
    # Distribution
    "BBB_Martins": {
        "type": "classification",
        "unit": "probability",
        "group": "distribution",
        "display_name": "BBB Penetration",
        "tdc_name": "BBB_Martins",
        "tdc_class": "ADME",
    },
    "PPBR_AZ": {
        "type": "regression",
        "unit": "%",
        "group": "distribution",
        "display_name": "Plasma Protein Binding",
        "tdc_name": "PPBR_AZ",
        "tdc_class": "ADME",
    },
    "VDss_Lombardo": {
        "type": "regression",
        "unit": "log L/kg",
        "group": "distribution",
        "display_name": "Volume of Distribution",
        "tdc_name": "VDss_Lombardo",
        "tdc_class": "ADME",
    },
    # Metabolism
    "CYP2D6_Veith": {
        "type": "classification",
        "unit": "probability",
        "group": "metabolism",
        "display_name": "CYP2D6 Inhibition",
        "tdc_name": "CYP2D6_Veith",
        "tdc_class": "ADME",
    },
    "CYP3A4_Veith": {
        "type": "classification",
        "unit": "probability",
        "group": "metabolism",
        "display_name": "CYP3A4 Inhibition",
        "tdc_name": "CYP3A4_Veith",
        "tdc_class": "ADME",
    },
    "CYP2C9_Veith": {
        "type": "classification",
        "unit": "probability",
        "group": "metabolism",
        "display_name": "CYP2C9 Inhibition",
        "tdc_name": "CYP2C9_Veith",
        "tdc_class": "ADME",
    },
    "CYP2D6_Substrate_CarbonMangels": {
        "type": "classification",
        "unit": "probability",
        "group": "metabolism",
        "display_name": "CYP2D6 Substrate",
        "tdc_name": "CYP2D6_Substrate_CarbonMangels",
        "tdc_class": "ADME",
    },
    "CYP3A4_Substrate_CarbonMangels": {
        "type": "classification",
        "unit": "probability",
        "group": "metabolism",
        "display_name": "CYP3A4 Substrate",
        "tdc_name": "CYP3A4_Substrate_CarbonMangels",
        "tdc_class": "ADME",
    },
    # Excretion
    "Half_Life_Obach": {
        "type": "regression",
        "unit": "hours",
        "group": "excretion",
        "display_name": "Half-life",
        "tdc_name": "Half_Life_Obach",
        "tdc_class": "ADME",
    },
    "Clearance_Hepatocyte_AZ": {
        "type": "regression",
        "unit": "uL/min/10^6 cells",
        "group": "excretion",
        "display_name": "Clearance (Hepatocyte)",
        "tdc_name": "Clearance_Hepatocyte_AZ",
        "tdc_class": "ADME",
    },
    "Clearance_Microsome_AZ": {
        "type": "regression",
        "unit": "uL/min/mg",
        "group": "excretion",
        "display_name": "Clearance (Microsome)",
        "tdc_name": "Clearance_Microsome_AZ",
        "tdc_class": "ADME",
    },
    # Toxicity
    "hERG": {
        "type": "classification",
        "unit": "probability",
        "group": "toxicity",
        "display_name": "hERG Cardiotoxicity",
        "tdc_name": "hERG",
        "tdc_class": "Tox",
    },
    "AMES": {
        "type": "classification",
        "unit": "probability",
        "group": "toxicity",
        "display_name": "AMES Mutagenicity",
        "tdc_name": "AMES",
        "tdc_class": "Tox",
    },
    "DILI": {
        "type": "classification",
        "unit": "probability",
        "group": "toxicity",
        "display_name": "Drug-Induced Liver Injury",
        "tdc_name": "DILI",
        "tdc_class": "Tox",
    },
    "LD50_Zhu": {
        "type": "regression",
        "unit": "log(1/(mol/kg))",
        "group": "toxicity",
        "display_name": "LD50 Acute Toxicity",
        "tdc_name": "LD50_Zhu",
        "tdc_class": "Tox",
    },
}
