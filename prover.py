from imported.fft import fft
from imported.kzg_proofs import is_power_of_two, get_root_of_unity, list_to_reverse_bit_order, commit_to_poly,\
    compute_proof_multi

from setup import get_setup
from shared import MODULUS, Sample, get_coset_factor
from my_types import G1Point


def create_matrix(blobs: [[int]], N_locs: int) -> ([[Sample]], [G1Point]):
    """Create sharding matrix and commitments from blobs"""

    # ensure data integrity
    assert all(len(blob) == len(get_setup()[0]) for blob in blobs)
    assert is_power_of_two(N_locs)
    assert len(blobs[0]) % N_locs == 0

    N_cols = len(blobs[0]) // N_locs

    matrix = []
    commitments = []

    # For each blob, compute the corresponding polynomial. Then split the blob into samples and compute their multiproofs.
    for i, blob in enumerate(blobs):
        polynomial = fft(list_to_reverse_bit_order(blob), MODULUS, get_root_of_unity(len(blobs[0])), True)

        row = []
        for j in range(N_cols):
            vs = blob[N_locs*j:N_locs*(j+1)]
            proof = compute_proof_multi(polynomial, get_coset_factor(j, N_locs), N_locs, get_setup())
            row.append(Sample(i, j, vs, proof))

        commitment = commit_to_poly(polynomial, get_setup())

        matrix.append(row)
        commitments.append(commitment)

    return matrix, commitments
