from py_ecc import optimized_bls12_381 as b
from py_ecc.fields import optimized_bls12_381_FQ as FQ
from py_ecc.typing import Optimized_Point3D

from imported.kzg_proofs import get_root_of_unity, list_to_reverse_bit_order

from setup import get_setup


MODULUS = b.curve_order


class Sample:
    i: int
    j: int
    vs: [int]
    proof: Optimized_Point3D[FQ]

    def __init__(self, i: int, j: int, vs: [int], proof: Optimized_Point3D[FQ]) -> None:
        self.i = i
        self.j = j
        self.vs = vs
        self.proof = proof


def get_coset_factor(j: int, N_locs: int) -> int:
    order = len(get_setup()[0])
    root = get_root_of_unity(order)
    return pow(root, list_to_reverse_bit_order(range(order))[N_locs * j], MODULUS)
