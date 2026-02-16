import csv
import io
import json
import os
import time
import uuid
from datetime import datetime, timedelta, timezone

from flask import Flask, jsonify, request, send_file
from flask_cors import CORS
from flask_jwt_extended import (
    JWTManager,
    create_access_token,
    get_jwt_identity,
    jwt_required,
)

from config.settings import (
    DATABASE_URL,
    JWT_ACCESS_TOKEN_EXPIRES_HOURS,
    JWT_SECRET_KEY,
    PREDICTIONS_CSV,
)
from models import User, db
from services.chemprop_service import ChempropService
from services.openai_service import OpenAIService
from services.report_service import ReportService
from services.training_service import TrainingService

app = Flask(__name__)
CORS(app, origins=["http://localhost:3000"])

# Database config
app.config["SQLALCHEMY_DATABASE_URI"] = DATABASE_URL
app.config["SQLALCHEMY_TRACK_MODIFICATIONS"] = False

# JWT config
app.config["JWT_SECRET_KEY"] = JWT_SECRET_KEY
app.config["JWT_ACCESS_TOKEN_EXPIRES"] = timedelta(hours=JWT_ACCESS_TOKEN_EXPIRES_HOURS)

db.init_app(app)
jwt = JWTManager(app)

with app.app_context():
    db.create_all()

chemprop_svc = ChempropService()
openai_svc = OpenAIService()
report_svc = ReportService()
training_svc = TrainingService()

# In-memory cache for report downloads (keyed by request_id)
_results_cache = {}


# ──────────────────────────── Auth Endpoints ────────────────────────────


@app.route("/api/auth/register", methods=["POST"])
def register():
    body = request.get_json(silent=True) or {}

    email = (body.get("email") or "").strip().lower()
    password = body.get("password") or ""
    first_name = (body.get("firstName") or "").strip()
    last_name = (body.get("lastName") or "").strip()
    organization = (body.get("organization") or "").strip() or None

    if not email or not password or not first_name or not last_name:
        return jsonify({"error": "Email, password, first name, and last name are required"}), 400

    if len(password) < 6:
        return jsonify({"error": "Password must be at least 6 characters"}), 400

    if User.query.filter_by(email=email).first():
        return jsonify({"error": "An account with this email already exists"}), 409

    user = User(
        email=email,
        first_name=first_name,
        last_name=last_name,
        organization=organization,
    )
    user.set_password(password)
    db.session.add(user)
    db.session.commit()

    token = create_access_token(identity=str(user.id))
    return jsonify({"token": token, "user": user.to_dict()}), 201


@app.route("/api/auth/login", methods=["POST"])
def login():
    body = request.get_json(silent=True) or {}

    email = (body.get("email") or "").strip().lower()
    password = body.get("password") or ""

    if not email or not password:
        return jsonify({"error": "Email and password are required"}), 400

    user = User.query.filter_by(email=email).first()
    if not user or not user.check_password(password):
        return jsonify({"error": "Invalid email or password"}), 401

    token = create_access_token(identity=str(user.id))
    return jsonify({"token": token, "user": user.to_dict()})


@app.route("/api/auth/me", methods=["GET"])
@jwt_required()
def me():
    user_id = get_jwt_identity()
    user = db.session.get(User, int(user_id))
    if not user:
        return jsonify({"error": "User not found"}), 404
    return jsonify({"user": user.to_dict()})


# ──────────────────────────── Health ────────────────────────────


@app.route("/api/health", methods=["GET"])
def health():
    return jsonify({
        "status": "ok",
        "mock_mode": chemprop_svc._use_mock,
        "models_loaded": chemprop_svc.models_loaded,
        "timestamp": datetime.now(timezone.utc).isoformat(),
    })


# ──────────────────────────── Predictions ────────────────────────────


@app.route("/api/predict", methods=["POST"])
@jwt_required()
def predict():
    start = time.time()
    request_id = str(uuid.uuid4())
    body = request.get_json(silent=True) or {}

    smiles_list = body.get("smiles", [])
    prompt = body.get("prompt", "")

    if not smiles_list:
        return jsonify({"error": "No SMILES provided"}), 400

    # Validate all SMILES
    validation = chemprop_svc.validate_smiles(smiles_list)
    if validation["invalid"]:
        return jsonify({
            "error": "Invalid SMILES",
            "invalid": validation["invalid"],
        }), 400

    # Run ADMET predictions for each molecule
    all_predictions = {}
    for smi in smiles_list:
        all_predictions[smi] = chemprop_svc.predict(smi)

    # GPT interpretation
    gpt_result = openai_svc.interpret(all_predictions, prompt)

    latency = round(time.time() - start, 2)

    result = {
        "request_id": request_id,
        "predictions": all_predictions,
        "interpretation": gpt_result,
        "mock_mode": chemprop_svc._use_mock,
        "latency_seconds": latency,
    }

    # Cache for report download
    _results_cache[request_id] = {
        "predictions": all_predictions,
        "interpretation": gpt_result,
        "smiles_list": smiles_list,
    }

    # Log to CSV
    _log_prediction(request_id, smiles_list, all_predictions, gpt_result, latency)

    return jsonify(result)


@app.route("/api/predict/<request_id>/report", methods=["GET"])
@jwt_required()
def download_report(request_id):
    cached = _results_cache.get(request_id)
    if not cached:
        return jsonify({"error": "Results expired or not found"}), 404

    report_path = report_svc.generate(
        request_id,
        cached["predictions"],
        cached["interpretation"],
        cached["smiles_list"],
    )
    return send_file(
        report_path,
        as_attachment=True,
        download_name=f"analysis_report_{request_id[:8]}.docx",
        mimetype="application/vnd.openxmlformats-officedocument.wordprocessingml.document",
    )


@app.route("/api/predict/<request_id>/csv", methods=["GET"])
@jwt_required()
def download_csv(request_id):
    """Generate and return a CSV file with prediction results for a request."""
    cached = _results_cache.get(request_id)
    if not cached:
        return jsonify({"error": "Results expired or not found"}), 404

    output = io.StringIO()
    writer = csv.writer(output)
    writer.writerow([
        "smiles", "endpoint", "group", "display_name",
        "value", "unit", "std", "type",
    ])

    for smiles, endpoints in cached["predictions"].items():
        for endpoint_name, data in endpoints.items():
            writer.writerow([
                smiles,
                endpoint_name,
                data.get("group", ""),
                data.get("display_name", endpoint_name),
                data.get("value", ""),
                data.get("unit", ""),
                data.get("std", ""),
                data.get("type", ""),
            ])

    csv_bytes = output.getvalue().encode("utf-8")
    output.close()

    return send_file(
        io.BytesIO(csv_bytes),
        as_attachment=True,
        download_name=f"admet_predictions_{request_id[:8]}.csv",
        mimetype="text/csv",
    )


# ──────────────────────────── Training ────────────────────────────


@app.route("/api/train", methods=["POST"])
def train():
    if training_svc.is_running:
        return jsonify({"error": "Training already in progress"}), 409

    body = request.get_json(silent=True) or {}
    endpoints = body.get("endpoints")  # optional: list of specific endpoints
    epochs = body.get("epochs")  # optional: override default

    training_svc.start(endpoints=endpoints, epochs=epochs)
    return jsonify({"status": "training_started"})


@app.route("/api/models/status", methods=["GET"])
def models_status():
    return jsonify(training_svc.get_status())


# ──────────────────────────── Helpers ────────────────────────────


def _log_prediction(request_id, smiles_list, predictions, gpt_result, latency):
    os.makedirs(os.path.dirname(PREDICTIONS_CSV), exist_ok=True)
    file_exists = os.path.exists(PREDICTIONS_CSV)
    with open(PREDICTIONS_CSV, "a", newline="") as f:
        writer = csv.writer(f)
        if not file_exists:
            writer.writerow([
                "timestamp_iso", "request_id", "smiles",
                "predictions_json", "ensemble_std_json",
                "gpt_summary_json", "latency_seconds", "status", "error_message",
            ])
        writer.writerow([
            datetime.now(timezone.utc).isoformat(),
            request_id,
            json.dumps(smiles_list),
            json.dumps(predictions),
            json.dumps({
                smi: {k: v.get("std") for k, v in preds.items()}
                for smi, preds in predictions.items()
            }),
            json.dumps(gpt_result.get("executive_summary", "")),
            latency,
            "success",
            "",
        ])


if __name__ == "__main__":
    app.run(debug=False, port=5001)
