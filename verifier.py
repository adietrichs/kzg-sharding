from py_ecc.fields import optimized_bls12_381_FQ as FQ
from py_ecc.typing import Optimized_Point3D

from imported.kzg_proofs import list_to_reverse_bit_order, check_proof_multi

from setup import get_setup
from shared import Sample, get_coset_factor


def verify(sample: Sample, commitments: [Optimized_Point3D[FQ]]) -> bool:
    """Verify a single sample multiproof"""
    commitment = commitments[sample.i]
    coset_factor = get_coset_factor(sample.j, len(sample.vs))
    ys = list_to_reverse_bit_order(sample.vs)
    return check_proof_multi(commitment, sample.proof, coset_factor, ys, get_setup())
