import pytest

from vexo.api.server import rdkit


@pytest.fixture()
def client():
    return rdkit.app.test_client()


def test_smart_substructure(client):
    resp = client.get("/smiles/substructure/CCC,C")
    assert resp.data == b"true\n"


def test_smiles_png(client):
    resp = client.get("/smiles/png/CCC")
    assert len(resp.data) == 3094
