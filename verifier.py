import random

from py_ecc import optimized_bls12_381 as b
from py_ecc.fields import optimized_bls12_381_FQ as FQ
from py_ecc.typing import Optimized_Point3D

from imported.fft import fft
from imported.kzg_proofs import div, list_to_reverse_bit_order, get_root_of_unity, check_proof_multi
from imported.multicombs import lincomb

from setup import get_setup
from shared import MODULUS, Sample, get_coset_factor


def verify(sample: Sample, commitments: [Optimized_Point3D[FQ]]) -> bool:
    """Verify a single sample multiproof"""
    commitment = commitments[sample.i]
    coset_factor = get_coset_factor(sample.j, len(sample.vs))
    ys = list_to_reverse_bit_order(sample.vs)
    return check_proof_multi(commitment, sample.proof, coset_factor, ys, get_setup())


def vector_entrywise_addition(vec_a: [int], vec_b: [int]) -> [int]:
    """Compute the entrywise addition of two vectors"""
    assert len(vec_a) == len(vec_b)
    return [a + b % MODULUS for a, b in zip(vec_a, vec_b)]


def verify_aggregated(samples: [Sample], commitments: [Optimized_Point3D[FQ]]) -> bool:
    """Verify multiple sample multiproofs at once"""
    if len(samples) == 0:
        return True

    N_rows = len(commitments)
    N_locs = len(samples[0].vs)
    N_cols = len(get_setup()[0]) // N_locs

    r = random.randint(1, MODULUS - 1)
    powers_of_r = [pow(r, k + 1, MODULUS) for k in range(len(samples))]

    proofs = [sample.proof for sample in samples]
    proof_lincomb = lincomb(proofs, powers_of_r, b.add, b.Z1)
    power_of_s = get_setup()[1][N_locs]

    commitment_weights = {i: 0 for i in range(N_rows)}
    for k, sample in enumerate(samples):
        commitment_weights[sample.i] += powers_of_r[k]
    used_commitments = [commitments[i] for i in range(N_rows) if commitment_weights[i] > 0]
    used_commitments_weights = [commitment_weights[i] for i in range(N_rows) if commitment_weights[i] > 0]
    g1_sum = lincomb(used_commitments, used_commitments_weights, b.add, b.Z1)

    aggregated_column_data = {j: [0] * N_locs for j in range(N_cols)}
    for sample, power_of_r in zip(samples, powers_of_r):
        scaled_data = [v * power_of_r for v in sample.vs] # scale the data points
        aggregated_column_sample[sample.j] = vector_entrywise_addition(aggregated_column_sample[sample.j], scaled_data)

    # We iterate over each sample and aggregate all the interpolation polynomials into `aggregated_interpolation_polynomial`
    aggregated_interpolation_polynomial = [0] * N_locs
    root_of_unity = get_root_of_unity(N_locs)
    for j in range(N_cols):
        if aggregated_column_data[j] == [0] * N_locs:
            continue
        coset_factor = get_coset_factor(j, N_locs)
        interpolation_polynomial = fft(list_to_reverse_bit_order(aggregated_column_data[j]), MODULUS, root_of_unity, True)
        interpolation_polynomial = [div(c, pow(coset_factor, i, MODULUS)) for i, c in enumerate(interpolation_polynomial)]
        # Update the aggregated interpolation polynomial
        aggregated_interpolation_polynomial = vector_entrywise_addition(aggregated_interpolation_polynomial, interpolation_polynomial)

    # Commit to the aggregated interpolation polynomial by evaluating it at `s` using the CRS
    evaluation = lincomb(
        get_setup()[0][:len(aggregated_interpolation_polynomial)], aggregated_interpolation_polynomial, b.add, b.Z1)
    g1_sum = b.add(g1_sum, b.neg(evaluation))

    weights = [pow(get_coset_factor(sample.j, N_locs), N_locs, MODULUS) for sample in samples]
    weighted_powers_of_r = [power_of_r * weight for power_of_r, weight in zip(powers_of_r, weights)]
    weighted_proof_lincomb = lincomb(proofs, weighted_powers_of_r, b.add, b.Z1)
    g1_sum = b.add(g1_sum, weighted_proof_lincomb)

    pairing_check = b.pairing(b.G2, b.neg(g1_sum), False)
    pairing_check *= b.pairing(power_of_s, proof_lincomb, False)
    pairing = b.final_exponentiate(pairing_check)
    return pairing == b.FQ12.one()
