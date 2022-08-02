import random

from py_ecc import optimized_bls12_381 as b

from imported.fft import fft
from imported.kzg_proofs import div, list_to_reverse_bit_order, get_root_of_unity
from imported.multicombs import lincomb

from setup import get_setup
from shared import MODULUS, Sample, get_coset_factor
from my_types import G1Point
from verifier import get_aggregated_pairings


def __detect_aggregated(samples: [Sample], commitments: [G1Point], begin: int, end: int) -> [int]:
    """ Detect if any samples/proofs in range [begin, end) of the row are corrupted"""
    # Obtain the pairings
    pairing_left, pairing_right = get_aggregated_pairings(samples[begin:end], commitments, begin + 1)

    # Do the pairing check
    pairing_check = pairing_left * pairing_right
    pairing = b.final_exponentiate(pairing_check)
    if pairing == b.FQ12.one():
        return []
    elif end == begin + 1:
        return [begin]
    
    # Split the samples and do binary search
    mid = (begin + end) // 2
    corrupted_list = __detect_aggregated(samples, commitments, begin, mid)
    # TODO: The correctness of another half can be obtained by
    #       pairing_left / pairing_left_of_left_samples == pairing_right / pairing_right_of_left_samples
    corrupted_list += __detect_aggregated(samples, commitments, mid, end)
    return corrupted_list


def detect_aggregated(samples: [Sample], commitments: [G1Point]) -> [int]:
    """
    Detect if any samples/proofs are corrupted using binary search.
    """
    return __detect_aggregated(samples, commitments, 0, len(samples))