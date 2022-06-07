from dataclasses import dataclass

from py_ecc import optimized_bls12_381 as b

from imported.kzg_proofs import get_root_of_unity, list_to_reverse_bit_order

from my_types import G1Point
from setup import get_setup


MODULUS = b.curve_order

@dataclass
class Sample():
    i: int
    j: int
    vs: [int]
    proof: G1Point

def get_coset_factor(j: int, N_locs: int) -> int:
    order = len(get_setup()[0])
    root = get_root_of_unity(order)
    return pow(root, list_to_reverse_bit_order(range(order))[N_locs * j], MODULUS)
