# taken from research/kzg_data_availability

from py_ecc import optimized_bls12_381 as b

from imported.fft import fft
from imported.multicombs import lincomb


PRIMITIVE_ROOT = 5
MODULUS = b.curve_order


def is_power_of_two(x):
    return x > 0 and x & (x-1) == 0


def get_root_of_unity(order):
    """
    Returns a root of unity of order "order"
    """
    assert (MODULUS - 1) % order == 0
    return pow(PRIMITIVE_ROOT, (MODULUS - 1) // order, MODULUS)


def inv(a):
    """
    Modular inverse using eGCD algorithm
    """
    if a == 0:
        return 0
    lm, hm = 1, 0
    low, high = a % MODULUS, MODULUS
    while low > 1:
        r = high // low
        nm, new = hm - lm * r, high - low * r
        lm, low, hm, high = nm, new, lm, low
    return lm % MODULUS


def div(x, y):
    return x * inv(y) % MODULUS


def div_polys(a, b):
    """
    Long polynomial difivion for two polynomials in coefficient form
    """
    a = [x for x in a]
    o = []
    apos = len(a) - 1
    bpos = len(b) - 1
    diff = apos - bpos
    while diff >= 0:
        quot = div(a[apos], b[bpos])
        o.insert(0, quot)
        for i in range(bpos, -1, -1):
            a[diff + i] -= b[i] * quot
        apos -= 1
        diff -= 1
    return [x % MODULUS for x in o]


def reverse_bit_order(n, order):
    """
    Reverse the bit order of an integer n
    """
    assert is_power_of_two(order)
    # Convert n to binary with the same number of bits as "order" - 1, then reverse its bit order
    return int(('{:0' + str(order.bit_length() - 1) + 'b}').format(n)[::-1], 2)


def list_to_reverse_bit_order(l):
    """
    Convert a list between normal and reverse bit order. This operation is idempotent.
    """
    return [l[reverse_bit_order(i, len(l))] for i in range(len(l))]


def commit_to_poly(polynomial, setup):
    """
    Kate commitment to polynomial in coefficient form
    """
    return lincomb(setup[0][:len(polynomial)], polynomial, b.add, b.Z1)


def compute_proof_multi(polynomial, x, n, setup):
    """
    Compute Kate proof for polynomial in coefficient form at positions x * w^y where w is
    an n-th root of unity (this is the proof for one data availability sample, which consists
    of several polynomial evaluations)
    """
    quotient_polynomial = div_polys(polynomial, [-pow(x, n, MODULUS)] + [0] * (n - 1) + [1])
    return lincomb(setup[0][:len(quotient_polynomial)], quotient_polynomial, b.add, b.Z1)


def check_proof_multi(commitment, proof, x, ys, setup):
    """
    Check a proof for a Kate commitment for an evaluation f(x w^i) = y_i
    """
    n = len(ys)
    root_of_unity = get_root_of_unity(n)

    # Interpolate at a coset. Note because it is a coset, not the subgroup, we have to multiply the
    # polynomial coefficients by x^i
    interpolation_polynomial = fft(ys, MODULUS, root_of_unity, True)
    interpolation_polynomial = [div(c, pow(x, i, MODULUS)) for i, c in enumerate(interpolation_polynomial)]

    # Verify the pairing equation
    #
    # e([commitment - interpolation_polynomial(s)], [1]) = e([proof],  [s^n - x^n])
    #    equivalent to
    # e([commitment - interpolation_polynomial]^(-1), [1]) * e([proof],  [s^n - x^n]) = 1_T
    #

    xn_minus_yn = b.add(setup[1][n], b.multiply(b.neg(b.G2), pow(x, n, MODULUS)))
    commitment_minus_interpolation = b.add(commitment, b.neg(lincomb(
        setup[0][:len(interpolation_polynomial)], interpolation_polynomial, b.add, b.Z1)))
    pairing_check = b.pairing(b.G2, b.neg(commitment_minus_interpolation), False)
    pairing_check *= b.pairing(xn_minus_yn, proof, False)
    pairing = b.final_exponentiate(pairing_check)
    return pairing == b.FQ12.one()
