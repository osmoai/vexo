from rdkit import Chem

from vexo.common import chemistry


def test_smiles_png():
    smiles = "C"
    assert chemistry.smiles_png(smiles) is not None


def test_canonicalize_smiles():
    smiles = "C1=CC=CC=C1"
    mol = Chem.MolFromSmiles(smiles)
    assert chemistry.canonicalize_smiles(mol) == "c1ccccc1"
