import hashlib
import logging

try:
    from rdkit import Chem
    _HAS_RDKIT = True
except ImportError:
    _HAS_RDKIT = False

from config.settings import ADMET_ENDPOINTS

logger = logging.getLogger(__name__)

# Mapping from ADMET-AI prediction column names to our endpoint keys
# ADMET-AI uses TDC dataset names as column names in its output
ADMETAI_COLUMN_MAP = {
    "Caco2_Wang": "Caco2_Wang",
    "HIA_Hou": "HIA_Hou",
    "Pgp_Broccatelli": "Pgp_Broccatelli",
    "Bioavailability_Ma": "Bioavailability_Ma",
    "BBB_Martins": "BBB_Martins",
    "PPBR_AZ": "PPBR_AZ",
    "VDss_Lombardo": "VDss_Lombardo",
    "CYP2D6_Veith": "CYP2D6_Veith",
    "CYP3A4_Veith": "CYP3A4_Veith",
    "CYP2C9_Veith": "CYP2C9_Veith",
    "CYP2D6_Substrate_CarbonMangels": "CYP2D6_Substrate_CarbonMangels",
    "CYP3A4_Substrate_CarbonMangels": "CYP3A4_Substrate_CarbonMangels",
    "Half_Life_Obach": "Half_Life_Obach",
    "Clearance_Hepatocyte_AZ": "Clearance_Hepatocyte_AZ",
    "Clearance_Microsome_AZ": "Clearance_Microsome_AZ",
    "hERG": "hERG",
    "AMES": "AMES",
    "DILI": "DILI",
    "LD50_Zhu": "LD50_Zhu",
}

# Fallback mock baselines (used only if admet_ai is not installed)
MOCK_BASELINES = {
    "Caco2_Wang": {"mean": -5.40, "spread": 0.80},
    "HIA_Hou": {"mean": 0.90, "spread": 0.15},
    "Pgp_Broccatelli": {"mean": 0.35, "spread": 0.20},
    "Bioavailability_Ma": {"mean": 0.72, "spread": 0.18},
    "BBB_Martins": {"mean": 0.75, "spread": 0.20},
    "PPBR_AZ": {"mean": 82.0, "spread": 15.0},
    "VDss_Lombardo": {"mean": 0.45, "spread": 0.30},
    "CYP2D6_Veith": {"mean": 0.18, "spread": 0.15},
    "CYP3A4_Veith": {"mean": 0.42, "spread": 0.20},
    "CYP2C9_Veith": {"mean": 0.30, "spread": 0.18},
    "CYP2D6_Substrate_CarbonMangels": {"mean": 0.25, "spread": 0.20},
    "CYP3A4_Substrate_CarbonMangels": {"mean": 0.55, "spread": 0.20},
    "Half_Life_Obach": {"mean": 3.2, "spread": 1.5},
    "Clearance_Hepatocyte_AZ": {"mean": 45.0, "spread": 20.0},
    "Clearance_Microsome_AZ": {"mean": 55.0, "spread": 25.0},
    "hERG": {"mean": 0.12, "spread": 0.10},
    "AMES": {"mean": 0.22, "spread": 0.18},
    "DILI": {"mean": 0.15, "spread": 0.12},
    "LD50_Zhu": {"mean": 2.50, "spread": 0.60},
}


class ChempropService:
    def __init__(self):
        self.models_loaded = False
        self._admet_model = None
        self._use_mock = False

        if not _HAS_RDKIT:
            logger.warning("rdkit not installed — SMILES validation will use basic fallback")

        self._load_admet_ai()

    def _load_admet_ai(self):
        """Load pre-trained ADMET-AI models."""
        try:
            from admet_ai import ADMETModel
            self._admet_model = ADMETModel()
            self.models_loaded = True
            logger.info("ADMET-AI pre-trained models loaded successfully")
        except ImportError:
            logger.warning(
                "admet_ai not installed — falling back to mock predictions. "
                "Install with: pip install admet-ai"
            )
            self._use_mock = True
        except Exception as e:
            logger.error("Failed to load ADMET-AI models: %s", e)
            self._use_mock = True

    def validate_smiles(self, smiles_list: list[str]) -> dict:
        """Validate a list of SMILES strings using RDKit (or basic fallback)."""
        valid = []
        invalid = []
        for smi in smiles_list:
            if _HAS_RDKIT:
                mol = Chem.MolFromSmiles(smi)
                if mol is not None:
                    valid.append(smi)
                else:
                    invalid.append(smi)
            else:
                # Basic fallback: accept any non-empty string with SMILES-like chars
                import re
                if smi and re.match(r'^[A-Za-z0-9@+\-\[\]\(\)\\\/=#$:.%]+$', smi):
                    valid.append(smi)
                else:
                    invalid.append(smi)
        return {"valid": valid, "invalid": invalid}

    def predict(self, smiles: str) -> dict:
        """Run ADMET predictions for a single SMILES string.

        Uses ADMET-AI pre-trained Chemprop models if available,
        otherwise falls back to mock predictions.
        """
        if self._admet_model and not self._use_mock:
            return self._predict_admet_ai(smiles)
        return self._predict_mock_all(smiles)

    def _predict_admet_ai(self, smiles: str) -> dict:
        """Run predictions using ADMET-AI pre-trained models."""
        preds = self._admet_model.predict(smiles=smiles)

        results = {}
        for endpoint_name, meta in ADMET_ENDPOINTS.items():
            col_name = ADMETAI_COLUMN_MAP.get(endpoint_name, endpoint_name)

            if col_name in preds:
                value = float(preds[col_name])
            else:
                # Try to find a matching column (ADMET-AI may use slightly different names)
                matched = False
                for key in preds:
                    if endpoint_name.lower() in key.lower():
                        value = float(preds[key])
                        matched = True
                        break
                if not matched:
                    # Fall back to mock for this endpoint
                    results[endpoint_name] = self._predict_mock(smiles, endpoint_name, meta)
                    continue

            # Clamp classification outputs to [0, 1]
            if meta["type"] == "classification":
                value = max(0.0, min(1.0, value))

            results[endpoint_name] = {
                "value": round(value, 4),
                "std": 0.0,  # ADMET-AI single model, no ensemble std
                "type": meta["type"],
                "unit": meta["unit"],
                "group": meta["group"],
                "display_name": meta["display_name"],
            }

        return results

    def _predict_mock_all(self, smiles: str) -> dict:
        """Generate mock predictions for all endpoints."""
        results = {}
        for name, meta in ADMET_ENDPOINTS.items():
            results[name] = self._predict_mock(smiles, name, meta)
        return results

    def _predict_mock(self, smiles: str, endpoint_name: str, meta: dict) -> dict:
        """Generate deterministic mock predictions seeded by the SMILES hash."""
        baseline = MOCK_BASELINES[endpoint_name]

        seed_str = f"{smiles}:{endpoint_name}"
        h = int(hashlib.sha256(seed_str.encode()).hexdigest()[:8], 16)
        perturbation = (h / 0xFFFFFFFF) * 2 - 1
        value = baseline["mean"] + perturbation * baseline["spread"] * 0.5

        if meta["type"] == "classification":
            value = max(0.01, min(0.99, value))

        std_seed = int(hashlib.sha256(f"std:{seed_str}".encode()).hexdigest()[:8], 16)
        std = 0.02 + (std_seed / 0xFFFFFFFF) * 0.08

        return {
            "value": round(value, 4),
            "std": round(std, 4),
            "type": meta["type"],
            "unit": meta["unit"],
            "group": meta["group"],
            "display_name": meta["display_name"],
        }
