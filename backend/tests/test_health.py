def test_health_returns_ok(client):
    resp = client.get("/api/health")
    assert resp.status_code == 200
    data = resp.get_json()
    assert data["status"] == "ok"
    assert "mock_mode" in data
    assert "timestamp" in data


def test_health_shows_mock_mode(client):
    resp = client.get("/api/health")
    data = resp.get_json()
    assert data["mock_mode"] is True
