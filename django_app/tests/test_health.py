def test_health_endpoint(api_client):
    r = api_client.get("/api/health/")
    assert r.status_code == 200
    assert r.json() == {"status": "ok", "service": "gsynth-api"}
