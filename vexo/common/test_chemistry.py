from vexo.common import chemistry


def test_smiles_png():
    smiles = "C"
    assert chemistry.smiles_png(smiles) is not None
