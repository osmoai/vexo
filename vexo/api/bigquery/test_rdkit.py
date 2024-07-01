import json

from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs

from vexo.api.bigquery import rdkit


class Request:
    def __init__(self, data):
        self._data = data

    def get_json(self):
        return self._data


def _make_request(calls):
    return Request({"calls": calls})


def test_iso_canonical_smiles():
    req = _make_request(
        [
            ["C"],
            ["CC1(C)[C@@H]2C[C@@H](O)[C@@](C)(O)[C@H]1C2"],
            ["Cc1cc(O)nc(C)n1"],
            ["CCO/C=N/c1ccccc1"],
        ]
    )
    resp = rdkit.iso_canonical_smiles(req)[0]
    resp = json.loads(resp)

    assert resp["replies"] == [c[0] for c in req.get_json()["calls"]]


def test_canonical_smiles():
    expected = ["C", "CC1(C)C2CC(O)C(C)(O)C1C2", "Cc1cc(O)nc(C)n1", "CCOC=Nc1ccccc1"]
    req = _make_request(
        [
            ["C"],
            ["CC1(C)[C@H]2C[C@@H]1[C@](C)(O)[C@H](O)C2"],
            ["Cc1cc(O)nc(C)n1"],
            ["CCO/C=N/c1ccccc1"],
        ]
    )
    resp = rdkit.canonical_smiles(req)[0]
    resp = json.loads(resp)

    assert resp["replies"] == expected


def test_inchi_inchikey():
    req = _make_request(
        [
            [
                "InChI=1S/C18H32O/c1-11-17(4,5)13-9-8-12-10-16(2,3)19-15(12)14(13)18(11,6)7/h11-15H,8-10H2,1-7H3"
            ],
            [
                "InChI=1S/C17H28O2/c1-10-15(2,3)11-8-9-12-14(13(11)16(10,4)5)19-17(6,7)18-12/h10,12,14H,8-9H2,1-7H3"
            ],
            [
                "InChI=1S/C14H22O/c1-9-13(2,3)10-7-6-8-11(15)12(10)14(9,4)5/h9H,6-8H2,1-5H3"
            ],
        ]
    )
    expected = [
        "HUYXPANWJVOYCI-UHFFFAOYSA-N",
        "QFQQCCDAGARSEN-UHFFFAOYSA-N",
        "MIZGSAALSYARKU-UHFFFAOYSA-N",
    ]
    resp = rdkit.inchi_inchikey(req)[0]
    resp = json.loads(resp)

    assert resp["replies"] == expected


def test_smiles_inchikey():
    req = _make_request(
        [
            ["CC1C(C)(C)C2CCC3CC(C)(C)OC3C2C1(C)C"],
            ["CC1C(C)(C)C2=C(C3OC(C)(C)OC3CC2)C1(C)C"],
            ["CC1C(C)(C)C2=C(C(=O)CCC2)C1(C)C"],
        ]
    )
    expected = [
        "HUYXPANWJVOYCI-UHFFFAOYSA-N",
        "QFQQCCDAGARSEN-UHFFFAOYSA-N",
        "MIZGSAALSYARKU-UHFFFAOYSA-N",
    ]
    resp = rdkit.smiles_inchikey(req)[0]
    resp = json.loads(resp)

    assert resp["replies"] == expected


def test_inchi_canonical_smiles():
    req = _make_request(
        [
            [
                "InChI=1S/C18H32O/c1-11-17(4,5)13-9-8-12-10-16(2,3)19-15(12)14(13)18(11,6)7/h11-15H,8-10H2,1-7H3"
            ],
            [
                "InChI=1S/C17H28O2/c1-10-15(2,3)11-8-9-12-14(13(11)16(10,4)5)19-17(6,7)18-12/h10,12,14H,8-9H2,1-7H3"
            ],
            [
                "InChI=1S/C14H22O/c1-9-13(2,3)10-7-6-8-11(15)12(10)14(9,4)5/h9H,6-8H2,1-5H3"
            ],
        ]
    )
    expected = [
        "CC1C(C)(C)C2CCC3CC(C)(C)OC3C2C1(C)C",
        "CC1C(C)(C)C2=C(C3OC(C)(C)OC3CC2)C1(C)C",
        "CC1C(C)(C)C2=C(C(=O)CCC2)C1(C)C",
    ]
    resp = rdkit.inchi_canonical_smiles(req)[0]
    resp = json.loads(resp)

    assert resp["replies"] == expected


def test_morgan_fingerprint():
    test_smiles = [
        ["C"],
        ["CC1(C)[C@H]2C[C@@H]1[C@](C)(O)[C@H](O)C2"],
        ["Cc1cc(O)nc(C)n1"],
        ["CCO/C=N/c1ccccc1"],
    ]
    test_mols = [[Chem.MolFromSmiles(smiles[0])] for smiles in test_smiles]

    req = _make_request(test_smiles)
    resp = rdkit.canonical_smiles(req)[0]
    resp = json.loads(resp)

    expected_fp = [
        AllChem.GetMorganFingerprintAsBitVect(
            mol[0], useChirality=True, radius=2, nBits=2048
        )
        for mol in test_mols
    ]
    expected = [DataStructs.BitVectToBinaryText(fp) for fp in expected_fp]

    resp = rdkit.morgan_fingerprint(req)[0]
    resp = json.loads(resp)

    test_results = [bytes.fromhex(result) for result in resp["replies"]]

    assert test_results == expected


def test_substructure_match():
    req = _make_request([["C", "CC1(C)[C@H]2C[C@@H](O)[C@@](C)(O)[C@H]1C2"]])
    resp = rdkit.substructure_match(req)[0]
    resp = json.loads(resp)
    assert resp["replies"] == [True]
