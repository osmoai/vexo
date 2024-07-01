import functions_framework

import vexo.api.server.rdkit as sheets_rdkit
import vexo.api.bigquery.rdkit as bq_rdkit


@functions_framework.http
def smart_substructure(request):
    query = request.path.lstrip("/")
    return sheets_rdkit.smart_substructure(query)


@functions_framework.http
def smiles_png(request):
    smiles = request.path.lstrip("/")
    return sheets_rdkit.smiles_png(smiles)


@functions_framework.http
def inchi_canonical_smiles(request):
    return bq_rdkit.inchi_canonical_smiles(request)


@functions_framework.http
def morgan_fingerprint(request):
    return bq_rdkit.morgan_ffn(request)


@functions_framework.http
def iso_canonical_smiles(request):
    return bq_rdkit.iso_canonical_smiles(request)


@functions_framework.http
def smiles_inchikey(request):
    return bq_rdkit.smiles_inchikey(request)


@functions_framework.http
def inchi_inchikey(request):
    return bq_rdkit.inchi_inchikey(request)


@functions_framework.http
def substructure_match(request):
    return bq_rdkit.substructure_match(request)


@functions_framework.http
def canonicalize_smiles(request):
    return bq_rdkit.canonicalize_smiles(request)
