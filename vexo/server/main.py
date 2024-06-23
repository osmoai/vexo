import argparse
import os

import flask
import waitress
from rdkit import Chem

import src.common.chemistry as chemistry

app = flask.Flask(__name__)


@app.route("/smiles/<path:smiles>")
def smiles(smiles: str):
    image_io = chemistry.smiles_png(smiles)
    return flask.send_file(image_io, mimetype="image/png")


@app.route("/smiles/substructure/<path:query>")
def smart_substructure(query: str):
    try:
        smiles, smart = query.split(",")
        mol = Chem.MolFromSmiles(smiles)
        fragment_mol = Chem.MolFromSmarts(smart)
        return flask.jsonify(mol.HasSubstructMatch(fragment_mol, useChirality=True))
    except Exception as e:
        return flask.abort(400, str(e))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Vexo")
    args = parser.parse_args()

    port = int(os.getenv("PORT", "8080"))
    waitress.serve(app, host="0.0.0.0", port=port, channel_timeout=3600)
