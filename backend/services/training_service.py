import json
import logging
import threading
from datetime import datetime, timezone

from config.settings import ADMET_ENDPOINTS, MODELS_DIR

logger = logging.getLogger(__name__)


class TrainingService:
    def __init__(self):
        self.is_running = False
        self._status = {
            "state": "idle",
            "current_endpoint": None,
            "completed_endpoints": [],
            "failed_endpoints": [],
            "started_at": None,
            "finished_at": None,
            "metrics": {},
        }

    def start(self, endpoints: list[str] | None = None, epochs: int | None = None):
        """Start training in a background thread."""
        self.is_running = True
        self._status = {
            "state": "running",
            "current_endpoint": None,
            "completed_endpoints": [],
            "failed_endpoints": [],
            "started_at": datetime.now(timezone.utc).isoformat(),
            "finished_at": None,
            "metrics": {},
        }
        thread = threading.Thread(
            target=self._run_training,
            args=(endpoints, epochs),
            daemon=True,
        )
        thread.start()

    def _run_training(self, endpoints: list[str] | None, epochs: int | None):
        from scripts.train_models import train_endpoint

        targets = endpoints or list(ADMET_ENDPOINTS.keys())

        for name in targets:
            if name not in ADMET_ENDPOINTS:
                continue

            self._status["current_endpoint"] = name
            logger.info("Training endpoint: %s", name)

            try:
                metrics = train_endpoint(name, ADMET_ENDPOINTS[name], epochs=epochs)
                self._status["completed_endpoints"].append(name)
                if metrics:
                    self._status["metrics"][name] = metrics
            except Exception as e:
                logger.error("Training failed for %s: %s", name, e)
                self._status["failed_endpoints"].append({
                    "name": name,
                    "error": str(e),
                })

        self.is_running = False
        self._status["state"] = "completed"
        self._status["current_endpoint"] = None
        self._status["finished_at"] = datetime.now(timezone.utc).isoformat()

    def get_status(self) -> dict:
        """Return current training status."""
        total = len(ADMET_ENDPOINTS)
        completed = len(self._status["completed_endpoints"])

        # Check which models have checkpoints on disk
        models_available = {}
        for name in ADMET_ENDPOINTS:
            model_dir = MODELS_DIR / name
            has_ckpts = model_dir.exists() and any(model_dir.glob("*.ckpt"))
            models_available[name] = has_ckpts

        # Load test metrics from disk if available
        saved_metrics = {}
        for name in ADMET_ENDPOINTS:
            metrics_path = MODELS_DIR / name / "test_metrics.json"
            if metrics_path.exists():
                with open(metrics_path) as f:
                    saved_metrics[name] = json.load(f)

        return {
            **self._status,
            "total_endpoints": total,
            "progress_percent": round(completed / total * 100) if total else 0,
            "models_available": models_available,
            "saved_metrics": saved_metrics,
        }
