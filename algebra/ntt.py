"""
The ntt module handles the Number Theoretic Transform (NTT) and its inverse in constant time.

Example usage:
    from ntt import find_primitive_root, bit_reverse_copy

    # Example parameters
    modulus = 17
    degree = 8  # we can infer this from the length of the input vector later!
    root_order = 2*degree

    # Precompute root and inv_root powers
    root = find_primitive_root(modulus=modulus, root_order=root_order)
    root_powers = [pow(base=root, exp=i, mod=modulus) for i in range(root_order)]
    brv_root_powers = bit_reverse_copy(val=root_powers)
    inv_root = pow(base=root, exp=modulus-2, mod=modulus)
    inv_root_powers = [pow(base=inv_root, exp=i, mod=modulus) for i in range(root_order)]
    brv_inv_root_powers = bit_reverse_copy(val=inv_root_powers)

    # Let's compute some NTTs!
    # Start with the multiplicative identity, the 1-polynomial
    one_poly = [1, 0, 0, 0, 0, 0, 0, 0]
    expected_one_poly_hat = [1, 1, 1, 1, 1, 1, 1, 1]  # manually computed
    observed_one_poly_hat = cooley_tukey_ntt(val=one_poly, modulus=modulus, root_order=root_order,
                                            bit_rev_root_powers=brv_root_powers)
    assert expected_one_poly_hat == observed_one_poly_hat

    # Next let's do the first monomial, the X-polynomial.
    x_poly = [0, 1, 0, 0, 0, 0, 0, 0]
    expected_x_poly_hat = [-8, 8, 8, -8, -2, 2, 2, -2]  # manually computed
    observed_x_poly_hat = cooley_tukey_ntt(val=x_poly, modulus=modulus, root_order=root_order,
                                           bit_rev_root_powers=brv_root_powers)
    assert expected_x_poly_hat == observed_x_poly_hat

    # Let's compute some inverse NTTs!
    # Let's check that we get back to the 1-polynomial from observed_one_poly_hat
    observed_one_poly_hat_hat = gentleman_sande(val=observed_one_poly_hat, modulus=modulus, root_order=root_order,
                                                bit_rev_inv_root_powers=brv_inv_root_powers)
    assert all((x-y) % modulus == 0 for x, y in zip(one_poly, one_poly_hat_hat))

"""
from typing import Dict, List

CACHED_IS_POW_TWO_GEQ_TWO: Dict[int, bool] = {}
_ERR_MUST_BE_INT: str = 'val must be an integer'


class InputValidationError(Exception):
    pass


def is_pow_two_geq_two(val: int) -> bool:
    """
    Check if a number is a power of two.
    :param val: Input integer
    :type val: int
    :return: True if val is a power of two, False otherwise
    :rtype: bool
    """
    if not isinstance(val, int):
        raise InputValidationError(_ERR_MUST_BE_INT)
    elif val not in CACHED_IS_POW_TWO_GEQ_TWO:
        CACHED_IS_POW_TWO_GEQ_TWO[val] = val >= 2 and val & (val - 1) == 0
    return CACHED_IS_POW_TWO_GEQ_TWO[val]


def bit_reverse_copy(val: list):
    """
    Permute indices by bit-reversal.
    :param val: Input list of anything
    :type val: list
    :return lav: Output list matches input list with indices bit-reversed
    :rtype: list
    """
    if not isinstance(val, list):
        raise TypeError("Input must be a list")
    elif not is_pow_two_geq_two(val=len(val)):
        raise ValueError("Input must be a list with power-of-two length.")
    n: int = len(val)
    k: int = n.bit_length() - 1
    bit_reversed_indices: List[int] = [int(bin(i)[2:].zfill(k)[::-1], 2) for i in range(n)]
    lav = [val[i] for i in bit_reversed_indices]  # bit reverse
    return lav


def cent(val: int, modulus: int, halfmod: int, logmod: int) -> int:
    """
    Centrally reduce a value modulo q in constant time. Output val satisfies -(q//2) <= val <= q//2.
    :param modulus: Input modulus
    :type modulus: int
    :param halfmod: Half the modulus
    :type halfmod: int
    :param logmod: bits to describe integers with this modulus
    :type logmod: int
    :param val: Input integer
    :type val: int
    :return z: such that z = val % q and -(q//2) <= z <= q//2
    :rtype: int
    """
    y: int = val % modulus
    intermediate_value: int = y - halfmod - 1
    z: int = y - (1 + (intermediate_value >> logmod)) * modulus
    return z


def cooley_tukey_ntt(val: List[int], modulus: int, root_order: int, bit_rev_root_powers: List[int], halfmod: int,
                     logmod: int) -> List[int]:
    """
    Input val, a list of n := len(val) integers in usual ordering, a modulus that is a prime such that
    (modulus-1) % (2*n) == 0, a root_order == 2*n, and a list of integers, root_powers, which are powers of a primitive
    root of unity with order root_order stored in bit-reversed order.  Output the NTT of val in bit-reversed order.

    In-place computation, iterative implementation of the Cooley-Tukey butterfly as defined in 'Speeding up the Number
    Theoretic Transform for Faster Ideal Lattice-Based Cryptography' by Longa and Naehrig
    (https://eprint.iacr.org/2016/504.pdf).

    :param val: Input list of integers
    :type val: List[int]
    :param modulus:
    :type modulus: int
    :param root_order:
    :type root_order: int
    :param bit_rev_root_powers:
    :type bit_rev_root_powers: List[int]
    :param halfmod: Half the modulus
    :type halfmod: int
    :param logmod: Bit length of modulus
    :type logmod: int
    :return val: Output list of integers
    :rtype: List[int]
    """
    u: int
    v: int
    n: int = len(val)
    t: int = n
    m: int = 1
    while m < len(val):
        t //= 2
        for i in range(m):
            j_one: int = 2 * i * t
            j_two: int = j_one + t - 1
            s: int = bit_rev_root_powers[m + i]
            for j in range(j_one, j_two + 1):
                u, v = val[j], val[j + t] * s
                val[j] = cent(val=u + v, modulus=modulus, halfmod=halfmod, logmod=logmod)
                val[j + t] = cent(val=u - v, modulus=modulus, halfmod=halfmod, logmod=logmod)
        m *= 2
    return val


def gentleman_sande_intt(val: List[int], modulus: int, root_order: int, bit_rev_inv_root_powers: List[int],
                         halfmod: int, logmod: int) -> List[int]:
    """
    Input val, a list of n := len(val) integers in usual ordering, a modulus that is a prime such that
    (modulus-1) % (2*n) == 0, a root_order == 2*n, and a list of integers, inv_root_powers, which are powers of the
    inverse of the primitive root of unity with order root_order used to compute the forward transform in bit-reversed
    order.  Output the INTT of val in standard order.

    In-place computation, iterative implementation of the Gentleman-Sande butterfly as in 'Speeding up the Number
    Theoretic Transform for Faster Ideal Lattice-Based Cryptography' by Longa and Naehrig
    (https://eprint.iacr.org/2016/504.pdf).

    :param val: Input list of integers
    :type val: List[int]
    :param modulus:
    :type modulus: int
    :param halfmod: Half the modulus
    :type halfmod: int
    :param logmod: Bit length of modulus
    :type logmod: int
    :param root_order:
    :type root_order: int
    :param bit_rev_inv_root_powers:
    :type bit_rev_inv_root_powers: List[int]
    :return val: Output list of integers
    :rtype: List[int]
    """
    u: int
    v: int
    n: int = len(val)
    n_inv: int = pow(n, modulus - 2, modulus)
    t: int = 1
    m: int = n
    while m > 1:
        j_one: int = 0
        h: int = m // 2
        for i in range(h):
            j_two: int = j_one + t - 1
            s: int = bit_rev_inv_root_powers[h + i]
            for j in range(j_one, j_two + 1):
                u, v = val[j], val[j + t]
                val[j] = cent(val=u + v, modulus=modulus, halfmod=halfmod, logmod=logmod)
                val[j + t] = cent(val=(u - v) * s, modulus=modulus, halfmod=halfmod, logmod=logmod)
            j_one += 2 * t
        t *= 2
        m = h
    for j in range(n):
        val[j] = cent(val=val[j] * n_inv, modulus=modulus, halfmod=halfmod, logmod=logmod)
    return val


def ntt_poly_mult(f: List[int], g: List[int], modulus: int, halfmod: int, logmod: int, degree: int, root_order: int,
                  root: int, inv_root: int) -> List[int]:
    """
    Multiply two coefficient representations of polynomials by first computing their NTTs
    and then their Hadamard product before inverting the NTT.
    """
    root_powers: List[int] = [pow(root, i, modulus) for i in range(degree)]
    brv_root_powers: List[int] = bit_reverse_copy(val=root_powers)
    inv_root_powers: List[int] = [pow(inv_root, i, modulus) for i in range(degree)]
    brv_inv_root_powers = bit_reverse_copy(val=inv_root_powers)
    cooley_tukey_ntt(val=f, modulus=modulus, root_order=root_order, bit_rev_root_powers=brv_root_powers,
                     halfmod=halfmod, logmod=logmod)
    cooley_tukey_ntt(val=g, modulus=modulus, root_order=root_order, bit_rev_root_powers=brv_root_powers,
                     halfmod=halfmod, logmod=logmod)
    fg: List[int] = [
        cent(val=x * y, modulus=modulus, halfmod=halfmod, logmod=logmod)
        for x, y in zip(f, g)
    ]
    gentleman_sande_intt(val=fg, modulus=modulus, root_order=root_order, bit_rev_inv_root_powers=brv_inv_root_powers,
                         halfmod=halfmod, logmod=logmod)
    gentleman_sande_intt(val=f, modulus=modulus, root_order=root_order, bit_rev_inv_root_powers=brv_inv_root_powers,
                         halfmod=halfmod, logmod=logmod)
    gentleman_sande_intt(val=g, modulus=modulus, root_order=root_order, bit_rev_inv_root_powers=brv_inv_root_powers,
                         halfmod=halfmod, logmod=logmod)
    return fg
