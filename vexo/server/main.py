import io

import flask
from rdkit import Chem

from PIL import Image
from rdkit import Chem
from rdkit.Chem import Draw

def smiles_png(request):
    # call rdkit to generate images
    image_io = io.BytesIO()
    smiles = request.path.lstrip('/')
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        img = Draw.MolToImage(mol)
    else:
        img = Image.new("RGB", (300, 300), (255, 255, 255))

    img.save(image_io, "PNG")
    image_io.seek(0)
    return flask.send_file(image_io, mimetype="image/png")


def smart_substructure(request):
    try:
        query = request.path.lstrip('/')
        smiles, smart = query.split(",")
        mol = Chem.MolFromSmiles(smiles)
        fragment_mol = Chem.MolFromSmarts(smart)
        return flask.jsonify(mol.HasSubstructMatch(fragment_mol, useChirality=True))
    except Exception as e:
        return flask.abort(400, str(e))
