"""Flask API for RDKit."""

import argparse
import os

import flask
from rdkit import Chem

from vexo.common import chemistry

app = flask.Flask(__name__)


@app.route("/smiles/png/<path:smiles>")
def smiles_png(smiles):
    """Render a PNG image of a molecule.
    Args:
        smiles: str, SMILES string.
    Returns:
        image/png: PNG image of the molecule.
    Examples:
        /smiles/png/CCC
        /smiles/png/C1=CC=CC=C1
    """
    image_io = chemistry.smiles_png(smiles)
    return flask.send_file(image_io, mimetype="image/png")


@app.route("/smiles/substructure/<path:query>")
def smart_substructure(query):
    """Check if a molecule contains a substructure.
    Args:
        query: str, comma-separated SMILES and SMARTS strings. Example: "CCC,C"
    Returns:
        bool: True if the molecule contains the substructure.
    Examples:
        /smiles/substructure/CCC,C
        /smiles/substructure/C1=CC=CC=C1,C1=CC=CC=C1
    """
    try:
        smiles, smart = query.split(",")
        mol = Chem.MolFromSmiles(smiles)
        fragment_mol = Chem.MolFromSmarts(smart)
        return flask.jsonify(mol.HasSubstructMatch(fragment_mol, useChirality=True))
    except Exception as e:
        return flask.abort(400, str(e))


if __name__ == "__main__":
    import waitress

    parser = argparse.ArgumentParser(description="Display")
    args = parser.parse_args()

    port = int(os.getenv("PORT", "8080"))
    waitress.serve(app, host="0.0.0.0", port=port, channel_timeout=3600)
