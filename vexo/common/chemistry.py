import functools
import io

from rdkit.Chem import Draw
from rdkit import Chem
from PIL import Image


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
