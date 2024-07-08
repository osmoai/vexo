import functools
import io

from PIL import Image
from rdkit import Chem
from rdkit.Chem import Draw


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
