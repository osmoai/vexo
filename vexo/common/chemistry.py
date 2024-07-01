import functools
import io

from PIL import Image
from rdkit import Chem
from rdkit.Chem import Draw


@functools.lru_cache(10_000)
def smiles_png(smiles: str):
    # call rdkit to generate images
    image_io = io.BytesIO()
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        img = Draw.MolToImage(mol)
    else:
        img = Image.new("RGB", (300, 300), (255, 255, 255))

    img.save(image_io, "PNG")
    image_io.seek(0)
    return image_io


def canonicalize_smiles(smiles: str, isomeric: bool = False):
    """Canonicalize SMILES string.
    Args:
        smiles: str, SMILES string.
        isomeric: bool, include stereochemistry.
    Returns:
        str, canonical SMILES string.
    Examples:
        canonical_smiles("C1=CC=CC=C1")
        canonical_smiles("C1=CC=CC=C1", isomeric=True)
    """
    mol = Chem.MolFromSmiles(smiles)
    return Chem.MolToSmiles(mol, isomericSmiles=isomeric, canonical=True)
