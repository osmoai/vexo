# Vexo

Make BigQuery and Google Sheets Chemistry Native

## Installation

Running `intall.sh` at the top directory will create cloud functions for all of the cheminformatic functions.  The caller
of the script must have access to create assign roles, create cloud functions, and create BigQuery datasets.  The default
project configured to the `gcloud` cli will be used to create the cloud and bigquery functions.

Currently the cloud functions targetting Google Sheets are not authenticated.  If you wish to authenticate them you will manually have to modify the `install.sh` script removing `--allow-unauthenticated` and adding a service account IAM binding like is done for the BigQuery functions.

## BigQuery

The `install.sh` script will create a dataset in BigQuery names `vexo`.  When it is complete there should be a function that matches each seen in `biquery_functions.csv`.

### Examples

#### Canonical Smiles

Canononicalize smiles based on standard rdkit canonicalization.

```
SELECT vexo.canonical_smiles("C1=CC=CC=C1")
```

#### Iso-Canonical Smiles

Same as above but keep stereochemistry in the returned smiles.

```
SELECT vexo.iso_canonical_smiles("CC1(C)[C@@H]2C[C@@H](O)[C@@](C)(O)[C@H]1C2")
```

### Substructure Search

Search for substructure match.  The first argument is the pattern and the second the smiles.  The pattern
can also be a SMARTS pattern.

```
SELECT vexo.substructure_match("c1ccccc1", "CCO/C=N/c1ccccc1")
```

### Morgan Fingerprints

Generate the morgan finterprints for a smiles.

```
SELECT vexo.morgan_fingerprint("CCO/C=N/c1ccccc1")
```


