import random
from unittest import TestCase

from py_ecc.fields import optimized_bls12_381_FQ as FQ
from py_ecc.typing import Optimized_Point3D

from prover import create_matrix
from setup import generate_setup
from shared import MODULUS, Sample
from verifier import verify


s = 1927409816240961209460912649124
N_rows = 2
N_cols = 2
N_locs = 16
generate_setup(s, N_cols * N_locs - 1)


class TestMatrix(TestCase):
    matrix: [[Sample]]
    commitments: [Optimized_Point3D[FQ]]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        blobs = [[random.randint(0, MODULUS - 1) for j in range(N_cols * N_locs)] for i in range(N_rows)]
        self.matrix, self.commitments = create_matrix(blobs, N_locs)

    def test_verify_sample_proofs(self):
        for row in self.matrix:
            for sample in row:
                self.assertTrue(verify(sample, self.commitments))
