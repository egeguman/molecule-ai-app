def test_predict_with_valid_smiles(client):
    resp = client.post("/api/predict", json={
        "smiles": ["CCO"],
        "prompt": "Analyze this molecule",
    })
    assert resp.status_code == 200
    data = resp.get_json()

    assert "request_id" in data
    assert "predictions" in data
    assert "interpretation" in data
    assert "CCO" in data["predictions"]
    assert data["mock_mode"] is True

    # Check ADMET structure
    preds = data["predictions"]["CCO"]
    assert "AMES" in preds
    assert "Caco2_Wang" in preds
    assert "value" in preds["AMES"]
    assert "std" in preds["AMES"]
    assert "type" in preds["AMES"]
    assert "group" in preds["AMES"]


def test_predict_multiple_smiles(client):
    resp = client.post("/api/predict", json={
        "smiles": ["CCO", "CC(=O)OC1=CC=CC=C1C(=O)O"],
        "prompt": "Compare these molecules",
    })
    assert resp.status_code == 200
    data = resp.get_json()
    assert len(data["predictions"]) == 2


def test_predict_invalid_smiles(client):
    resp = client.post("/api/predict", json={
        "smiles": ["NOT_A_VALID_SMILES"],
        "prompt": "test",
    })
    assert resp.status_code == 400
    data = resp.get_json()
    assert "invalid" in data


def test_predict_empty_smiles_list(client):
    resp = client.post("/api/predict", json={
        "smiles": [],
        "prompt": "",
    })
    assert resp.status_code == 400


def test_predict_no_body(client):
    resp = client.post("/api/predict")
    assert resp.status_code == 400


def test_predict_deterministic_mock(client):
    """Same SMILES should produce the same mock predictions."""
    resp1 = client.post("/api/predict", json={"smiles": ["CCO"], "prompt": ""})
    resp2 = client.post("/api/predict", json={"smiles": ["CCO"], "prompt": ""})

    preds1 = resp1.get_json()["predictions"]["CCO"]
    preds2 = resp2.get_json()["predictions"]["CCO"]

    for endpoint in preds1:
        assert preds1[endpoint]["value"] == preds2[endpoint]["value"]


def test_predict_interpretation_structure(client):
    resp = client.post("/api/predict", json={
        "smiles": ["CCO"],
        "prompt": "Analyze ADMET",
    })
    data = resp.get_json()
    interp = data["interpretation"]

    assert "executive_summary" in interp
    assert "high_risk_flags" in interp
    assert "recommended_next_steps" in interp
    assert "disclaimer" in interp


def test_models_status(client):
    resp = client.get("/api/models/status")
    assert resp.status_code == 200
    data = resp.get_json()
    assert "state" in data
    assert "total_endpoints" in data
    assert data["total_endpoints"] == 19
