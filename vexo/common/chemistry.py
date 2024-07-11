import functools
import io

from PIL import Image
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.MolStandardize import rdMolStandardize


@functools.cache
def _tautomer_enumeration():
    enumerator = rdMolStandardize.TautomerEnumerator()
    enumerator.SetRemoveSp3Stereo(False)
    enumerator.SetRemoveBondStereo(False)
    return enumerator


def smiles_png(smiles: str):
    """Render a PNG image of a molecule.
    Args:
        smiles: str, SMILES string.
    Returns:
        io.BytesIO: PNG image of the molecule.
    Examples:
        smiles_png("CCC")
        smiles_png("C1=CC=CC=C1")
    """
    image_io = io.BytesIO()
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        img = Draw.MolToImage(mol)
    else:
        img = Image.new("RGB", (300, 300), (255, 255, 255))

    img.save(image_io, "PNG")
    image_io.seek(0)
    return image_io


def canonicalize_smiles(mol: Chem.Mol, isomeric: bool = False) -> str:
    """Canonicalize SMILES string.
    Args:
        mol: Mol, RDKit molecule.
        isomeric: bool, include stereochemistry.
    Returns:
        str, canonical SMILES string.
    """
    enumerator = _tautomer_enumeration()
    mol = enumerator.Canonicalize(mol)
    return Chem.MolToSmiles(mol, isomericSmiles=isomeric, canonical=True)
