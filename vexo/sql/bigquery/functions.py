import functools

from google.cloud import bigquery


@functools.cache
def _client():
    return bigquery.Client()


def _cos_similarity():
    _client().query(
        """
        CREATE OR REPLACE FUNCTION `vexo.cos_similarity`(
            arr1 ARRAY<FLOAT64>,
            arr2 ARRAY<FLOAT64>
        ) AS (
            (
                SELECT
                    CASE
                    WHEN ARRAY_LENGTH(arr1) != ARRAY_LENGTH(arr2) THEN NULL
                    WHEN arr1 is NULL OR arr2 is NULL THEN NULL
                    WHEN ARRAY_LENGTH(arr1) < 1 THEN NULL
                    ELSE ROUND(vexo.dot(vexo.norm_vector(arr1), vexo.norm_vector(arr2)), 2)
                END
            )
        )
        """
    ).result()


def _cos_distance():
    _client().query(
        """
        CREATE OR REPLACE FUNCTION `vexo.cos_distance`(
            arr1 ARRAY<FLOAT64>,
            arr2 ARRAY<FLOAT64>
        ) AS (
            (SELECT 1. - vexo.cos_similarity(arr1, arr2))
        )
        """
    ).result()


def _dot():
    _client().query(
        """
        CREATE OR REPLACE FUNCTION `vexo.dot`(
            arr1 ARRAY<FLOAT64>,
            arr2 ARRAY<FLOAT64>
        ) AS (
            (
                SELECT SUM(arr1[OFFSET(idx)] * arr2[OFFSET(idx)])
                FROM UNNEST(GENERATE_ARRAY(0, ARRAY_LENGTH(arr1) - 1)) as idx
            )
        )
        """
    ).result()


def _l1_distance():
    _client().query(
        """
        CREATE OR REPLACE FUNCTION `vexo.l1_distance`(
            arr1 ARRAY<FLOAT64>,
            arr2 ARRAY<FLOAT64>
        ) AS (
            (
                SELECT SUM(ABS(arr1[OFFSET(idx)] - arr2[OFFSET(idx)]))
                FROM UNNEST(GENERATE_ARRAY(0, ARRAY_LENGTH(arr1) - 1)) as idx
            )
        )
        """
    ).result()


def _l2_distance():
    _client().query(
        """
        CREATE OR REPLACE FUNCTION `vexo.l2_distance`(
            arr1 ARRAY<FLOAT64>,
            arr2 ARRAY<FLOAT64>
        ) AS (
            (
                SELECT SQRT(SUM(POW((arr1[OFFSET(idx)] - arr2[OFFSET(idx)]), 2)))
                FROM UNNEST(GENERATE_ARRAY(0, ARRAY_LENGTH(arr1) - 1)) as idx
            )
        )
        """
    ).result()


def _norm_vector():
    _client().query(
        """
        CREATE OR REPLACE FUNCTION `vexo.norm_vector`(
            arr1 ARRAY<FLOAT64>
        ) AS (
            (SELECT ARRAY_AGG(a/vexo.vector_length(arr1)) FROM UNNEST(arr1) as a)
        )
        """
    ).result()


def _tanimoto_hex():
    _client().query(
        """
        CREATE OR REPLACE FUNCTION `vexo.tanimoto_hex`(
            arr1 STRING,
            arr2 STRING
        ) AS (
            (1 - (BIT_COUNT(FROM_HEX(arr1) & FROM_HEX(arr2)) / BIT_COUNT(FROM_HEX(arr1) | FROM_HEX(arr2))))
        )
        """
    ).result()


def _vector_length():
    _client().query(
        """
        CREATE OR REPLACE FUNCTION `vexo.vector_length`(
            arr ARRAY<FLOAT64>
        ) AS (
            (SELECT SQRT(SUM(POW(a,2))) FROM UNNEST(arr) as a)
        )
        """
    ).result()


def main():
    _vector_length()
    _norm_vector()
    _dot()
    _l1_distance()
    _l2_distance()
    _tanimoto_hex()
    _cos_similarity()
    _cos_distance()


if __name__ == "__main__":
    main()
