"""BigQuery Cloud Function around rdkit"""

import json

from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs


def _is_valid_smiles(smi: str) -> bool:
    """Check if SMILES is valid"""
    mol = Chem.MolFromSmiles(smi, sanitize=False)
    if mol is None:
        return False
    try:
        Chem.SanitizeMol(mol)
        return True
    except ValueError:
        return False


def _bq_fn(request, fn):
    """Common function for BigQuery Cloud Function
    Args:
        request: Flask request object
        fn: function to apply to each call
    Returns:
        json response
    """
    try:
        return_value = []
        request_json = request.get_json()
        calls = request_json["calls"]
        for call in calls:
            try:
                return_value.append(fn(call))
            except Exception:
                return_value.append("")

        return_json = json.dumps({"replies": return_value}), 200
        return return_json
    except Exception:
        return json.dumps({"errorMessage": "something unexpected in input"}), 400


def canonical_smiles(request):
    """Canonicalize SMILES in BigQuery"""

    def fn(call):
        smiles = call[0]
        if _is_valid_smiles(smiles):
            mol = Chem.MolFromSmiles(smiles)
            result = Chem.MolToSmiles(mol, isomericSmiles=False, canonical=True)
            return result
        else:
            return ""

    return _bq_fn(request, fn)


def inchi_inchikey(request):
    """Inchi to Inchikey"""

    def fn(call):
        return Chem.inchi.InchiToInchiKey(call[0])

    return _bq_fn(request, fn)


def smiles_inchikey(request):
    """Smiles to Inchikey"""

    def fn(call):
        smiles = call[0]
        mol = Chem.MolFromSmiles(smiles)
        return Chem.MolToInchiKey(mol)

    return _bq_fn(request, fn)


def iso_canonical_smiles(request):
    """Canonicalize SMILES in BigQuery keeping stereochemistry"""

    def fn(call):
        smiles = call[0]
        if _is_valid_smiles(smiles):
            mol = Chem.MolFromSmiles(smiles)
            result = Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
            return result
        else:
            return ""

    return _bq_fn(request, fn)


def morgan_fingerprint(request):
    """Converts SMILES to Morgan fingerprint in BigQuery"""

    def fn(call):
        smiles = call[0]
        mol = Chem.MolFromSmiles(smiles)
        fp = AllChem.GetMorganFingerprintAsBitVect(
            mol, useChirality=True, radius=2, nBits=2048
        )
        return DataStructs.BitVectToBinaryText(fp).hex()

    return _bq_fn(request, fn)


def inchi_canonical_smiles(request):
    """Convert inchi to canonical SMILES in BigQuery"""

    def fn(call):
        inchi = call[0]
        mol = Chem.MolFromInchi(inchi)
        return Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)

    return _bq_fn(request, fn)


def substructure_match(request):
    """Check if a fragment is a substructure of a molecule"""
    try:
        return_value = []
        request_json = request.get_json()
        calls = request_json["calls"]

        for call in calls:

            fragment_smiles = call[0]
            smiles = call[1]

            try:
                mol = Chem.MolFromSmiles(smiles)
                fragment_mol = Chem.MolFromSmarts(fragment_smiles)
                return_value.append(
                    mol.HasSubstructMatch(fragment_mol, useChirality=True)
                )

            except Exception:
                return_value.append(False)

        return_json = json.dumps({"replies": return_value}), 200
        return return_json
    except Exception:
        return json.dumps({"errorMessage": "something unexpected in input"}), 400
