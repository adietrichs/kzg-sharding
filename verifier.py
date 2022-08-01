import random

from py_ecc import optimized_bls12_381 as b

from imported.fft import fft
from imported.kzg_proofs import div, list_to_reverse_bit_order, get_root_of_unity, check_proof_multi
from imported.multicombs import lincomb

from setup import get_setup
from shared import MODULUS, Sample, get_coset_factor
from my_types import G1Point


def verify(sample: Sample, commitments: [G1Point]) -> bool:
    """Verify a single sample multiproof"""
    commitment = commitments[sample.i]
    coset_factor = get_coset_factor(sample.j, len(sample.vs))
    ys = list_to_reverse_bit_order(sample.vs)
    return check_proof_multi(commitment, sample.proof, coset_factor, ys, get_setup())


def vector_entrywise_addition(vec_a: [int], vec_b: [int]) -> [int]:
    """Compute the entrywise addition of two vectors"""
    assert len(vec_a) == len(vec_b)
    return [a + b % MODULUS for a, b in zip(vec_a, vec_b)]


def get_aggregated_pairings(samples: [Sample], commitments: [G1Point], rpower_base: int) -> (b.FQ12, b.FQ12):
    """Verify multiple sample multiproofs at once"""
    assert len(samples) > 0

    N_rows = len(commitments)  # Number of rows of the blob matrix
    N_locs = len(samples[0].vs)  # Number of data points in each sample
    N_cols = len(get_setup()[0]) // N_locs  # Number of columns of the blob matrix

    # Derive random factors for the random linear combination
    r = random.randint(1, MODULUS - 1)
    powers_of_r = [pow(r, k, MODULUS) for k in range(rpower_base, rpower_base + len(samples))]

    # Here is a (simplified version) of the verification formula:
    #     e(∑ₖ rᵏ πₖ, [s¹⁶]₂) = e(∑ᵢ(∑ rᵏ)Cᵢ − [I(s)]₁ + ∑ₖ rᵏ hⱼ¹⁶πₖ, [1]₂)
    #        where I(x) = ∑ₖ rᵏ Iₖ(x) = ∑ⱼ∑ₗ(∑ rᵏ νₖₗ)Lⱼₗ(x)
    # In the following code we will be computing the verification formula from left-to-right:

    # Step 1) Compute random linear combination of the proofs
    proofs = [sample.proof for sample in samples]
    proof_lincomb = lincomb(proofs, powers_of_r, b.add, b.Z1)

    # Step 2) Get [s^16]
    power_of_s = get_setup()[1][N_locs]

    # Step 3) Compute sum of the commitments
    # First compute sum of the random lincomb factors for each commitment
    commitment_weights = [0] * N_rows
    for k, sample in enumerate(samples):
        commitment_weights[sample.i] += powers_of_r[k]
    # Find the commitments that correspond to the received samples
    used_commitments = [commitments[i] for i in range(N_rows) if commitment_weights[i] > 0]
    used_commitments_weights = [commitment_weights[i] for i in range(N_rows) if commitment_weights[i] > 0]
    # Compute commitment sum
    final_g1_sum = lincomb(used_commitments, used_commitments_weights, b.add, b.Z1)

    # Step 4) Compute sum of the interpolation polynomials
    # To do this, we perform the following logic:
    # a) For each column, aggregate all data into an aggregated column sample (stored in `aggregated_column_samples`)
    # b) For each column, compute interpolation polynomial (in coefficient form) that corresponds to the aggregated column sample
    # c) Finally, add all column interpolation polynomials together to get the final aggregated interpolation polynomial

    # We use `aggregated_column_samples` to store the scaled data points of each column
    aggregated_column_samples = {j: [0] * N_locs for j in range(N_cols)}
    # Iterate over each sample, scale its data points and update `aggregated_column_data`
    for sample, power_of_r in zip(samples, powers_of_r):
        scaled_data = [v * power_of_r for v in sample.vs] # scale the data points
        aggregated_column_samples[sample.j] = vector_entrywise_addition(aggregated_column_samples[sample.j], scaled_data)

    # The final aggregated interpolation polynomial
    aggregated_interpolation_poly = [0] * N_locs
    # Iterate over each column and compute the column interpolation polynomial
    root_of_unity = get_root_of_unity(N_locs) # generator for a small roots of unity subgroup used for interpolations
    for j in range(N_cols):
        if aggregated_column_samples[j] == [0] * N_locs: # skip this column if we have no samples from it
            continue
        coset_factor = get_coset_factor(j, N_locs) # get coset factor for this column

        # Get interpolation polynomial for this column. To do so we first do an IDFT over the roots of unity and then we
        # scale by the coset factor. We can't do an IDFT directly over the coset because it's not a subgroup.
        column_interpolation_poly = fft(list_to_reverse_bit_order(aggregated_column_samples[j]), MODULUS, root_of_unity, True)
        column_interpolation_poly = [div(c, pow(coset_factor, i, MODULUS)) for i, c in enumerate(column_interpolation_poly)]

        # While we are at it, update the final aggregated interpolation polynomial
        aggregated_interpolation_poly = vector_entrywise_addition(aggregated_interpolation_poly, column_interpolation_poly)

    # Commit to the final aggregated interpolation polynomial
    evaluation = lincomb(
        get_setup()[0][:len(aggregated_interpolation_poly)], aggregated_interpolation_poly, b.add, b.Z1)
    final_g1_sum = b.add(final_g1_sum, b.neg(evaluation))

    # Step 5) Compute sum of the proofs scaled by the coset factors
    weights = [pow(get_coset_factor(sample.j, N_locs), N_locs, MODULUS) for sample in samples]
    weighted_powers_of_r = [power_of_r * weight for power_of_r, weight in zip(powers_of_r, weights)]
    weighted_proof_lincomb = lincomb(proofs, weighted_powers_of_r, b.add, b.Z1)
    final_g1_sum = b.add(final_g1_sum, weighted_proof_lincomb)

    # Step 6) Calculate the final pairings
    pairing_right = b.pairing(b.G2, b.neg(final_g1_sum), False)
    pairing_left = b.pairing(power_of_s, proof_lincomb, False)
    return pairing_left, pairing_right


def verify_aggregated(samples: [Sample], commitments: [G1Point]) -> bool:
    """Verify multiple sample multiproofs at once"""
    if len(samples) == 0:
        return True

    # Obtain the pairings
    pairing_left, pairing_right = get_aggregated_pairings(samples, commitments, 1)

    # Do the final pairing check
    pairing_check = pairing_left * pairing_right
    pairing = b.final_exponentiate(pairing_check)
    return pairing == b.FQ12.one()
