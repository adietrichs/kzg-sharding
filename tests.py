import random
from unittest import TestCase

from py_ecc.fields import optimized_bls12_381_FQ as FQ

from prover import create_matrix
from setup import generate_setup
from shared import MODULUS, Sample
from verifier import verify, verify_aggregated

from my_types import G1Point

s = 1927409816240961209460912649124
N_rows = 4
N_cols = 4
N_locs = 16
generate_setup(s, N_cols * N_locs - 1)


class TestMatrix(TestCase):
    matrix: [[Sample]]
    commitments: [G1Point]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        blobs = [[random.randint(0, MODULUS - 1) for j in range(N_cols * N_locs)] for i in range(N_rows)]
        self.matrix, self.commitments = create_matrix(blobs, N_locs)

    def test_verify_sample_proofs(self):
        for i in range(1, 3):
            for j in range(2, 4):
                self.assertTrue(verify(self.matrix[i][j], self.commitments))

    def test_verify_aggregated_sample_proofs(self):
        samples = [self.matrix[0][3], self.matrix[2][0], self.matrix[2][2], self.matrix[3][2]]
        self.assertTrue(verify_aggregated(samples, self.commitments))

if __name__ == '__main__':
    import unittest
    unittest.main()
