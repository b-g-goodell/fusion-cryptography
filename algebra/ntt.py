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
from math import ceil
from typing import Dict, Tuple, List, Union
from algebra.errors import MUST_BE_INT_ERR, MUST_BE_LIST_ERR, MUST_BE_LIST_W_POW_2_LEN_ERR, MUST_BE_POS_INT_ERR, MUST_BE_INT_GEQ_3_ERR, MUST_BE_HALF_FLOORED_ERR, MUST_BE_BIT_LEN_ERR, MUST_HAVE_PROU_ERR, NO_PROU_FOUND_ERR, MUST_BE_BOOL_ERR, MUST_BE_ODD_PRIME_ERR, MUST_BE_CORRECT_ROOT_ERR, MUST_BE_CORRECT_INVERSE_ROOT_ERR, DEGREE_MISMATCH_ERR, MUST_CONSTRUCT_BRV_POWERS_CORRECTLY_ERR


# Caching dictionaries
CACHED_PRIMITIVE_ROOTS: Dict[Tuple[int, int], int] = {}
CACHED_IS_ODD_PRIME: Dict[int, bool] = {}
CACHED_HAS_PRIMITIVE_ROOT_OF_UNITY: Dict[Tuple[int, int], bool] = {}
CACHED_IS_POW_TWO_GEQ_TWO: Dict[int, bool] = {}
CACHED_IS_ROOT_OF_UNITY: Dict[Tuple[int, int, int], bool] = {}
CACHED_IS_PRIMITIVE_ROOT_OF_UNITY: Dict[Tuple[int, int, int], bool] = {}
CACHED_FIND_PRIMITIVE_ROOT: Dict[Tuple[int, int], Union[int, None]] = {}
CACHED_REVERSED_INDICES = {}


# Functions
def is_odd_prime(val: int) -> bool:
    """
    Check if x is an odd prime number.
    """
    if not isinstance(val, int):
        raise TypeError(MUST_BE_INT_ERR)
    elif val not in CACHED_IS_ODD_PRIME:
        CACHED_IS_ODD_PRIME[val] = val >= 3 and all(val % i != 0 for i in range(2, ceil(val**0.5) + 1))
    return CACHED_IS_ODD_PRIME[val]


def has_prou(mod: int, deg: int) -> bool:
    """
    Check if a modulus has a primitive root of unity of a given order.
    """
    if not isinstance(mod, int) or not isinstance(deg, int):
        raise TypeError(MUST_BE_INT_ERR)
    elif (mod, deg) not in CACHED_HAS_PRIMITIVE_ROOT_OF_UNITY:
        CACHED_HAS_PRIMITIVE_ROOT_OF_UNITY[(mod, deg)] = mod >= 3 and deg >= 1 and (mod - 1) % (2*deg) == 0
    return CACHED_HAS_PRIMITIVE_ROOT_OF_UNITY[(mod, deg)]


def is_pow_two_geq_two(val: int) -> bool:
    """
    Check if a number is a power of two.
    """
    if not isinstance(val, int):
        raise TypeError(MUST_BE_INT_ERR)
    elif val not in CACHED_IS_POW_TWO_GEQ_TWO:
        CACHED_IS_POW_TWO_GEQ_TWO[val] = val >= 2 and val & (val - 1) == 0
    return CACHED_IS_POW_TWO_GEQ_TWO[val]


def bit_reverse_copy(val: list):
    """
    Permute indices by bit-reversal.
    """
    if not isinstance(val, list):
        raise TypeError(MUST_BE_LIST_ERR)
    elif not is_pow_two_geq_two(val=len(val)):
        raise ValueError(MUST_BE_LIST_W_POW_2_LEN_ERR)
    elif len(val) not in CACHED_REVERSED_INDICES:
        k: int = len(val).bit_length() - 1
        CACHED_REVERSED_INDICES[len(val)] = [int(bin(i)[2:].zfill(k)[::-1], 2) for i in range(len(val))]
    return [val[i] for i in CACHED_REVERSED_INDICES[len(val)]]


def check_mod_halfmod_logmod(mod: int, halfmod: int, logmod: int):
    if any(not isinstance(x, int) for x in [mod, halfmod, logmod]):
        raise TypeError(MUST_BE_INT_ERR)
    elif mod < 1 or halfmod < 1 or logmod < 1:
        raise ValueError(MUST_BE_POS_INT_ERR)
    elif mod < 3:
        raise ValueError(MUST_BE_INT_GEQ_3_ERR)
    elif halfmod != mod//2:
        raise ValueError(MUST_BE_HALF_FLOORED_ERR)
    elif logmod != mod.bit_length():
        raise ValueError(MUST_BE_BIT_LEN_ERR)
    pass


def is_rou(root: int, mod: int, deg: int) -> bool:
    """
    Check if val is a root of unity of order root_order modulo modulus.
    """
    if not isinstance(root, int) or not isinstance(mod, int) or not isinstance(deg, int):
        raise TypeError(MUST_BE_INT_ERR)
    elif (root, mod, deg) not in CACHED_IS_ROOT_OF_UNITY:
        CACHED_IS_ROOT_OF_UNITY[(root, mod, deg)] = mod >= 3 and deg >= 1 and pow(base=root, exp=2*deg, mod=mod) == 1
    return CACHED_IS_ROOT_OF_UNITY[(root, mod, deg)]


def is_prou(root: int, mod: int, deg: int) -> bool:
    """
    Check if val is a primitive root of order root_order modulo modulus.
    """
    if not isinstance(root, int) or not isinstance(mod, int) or not isinstance(deg, int):
        raise TypeError(MUST_BE_INT_ERR)
    elif (root, mod, deg) not in CACHED_IS_PRIMITIVE_ROOT_OF_UNITY:
        is_prim: bool = is_rou(root=root, mod=mod, deg=deg)
        is_prim = is_prim and all(pow(base=root, exp=i, mod=mod) != 1 for i in range(1, 2*deg))
        CACHED_IS_PRIMITIVE_ROOT_OF_UNITY[(root, mod, deg)] = is_prim
    return CACHED_IS_PRIMITIVE_ROOT_OF_UNITY[(root, mod, deg)]


def find_prou(mod: int, deg: int) -> int:
    """
    Find a primitive root of order root_order modulo modulus. Naive loop that first checks 2, then 3, then 4...
    """
    if not isinstance(mod, int) or not isinstance(deg, int):
        raise TypeError(MUST_BE_INT_ERR)
    elif not has_prou(mod=mod, deg=deg):
        raise ValueError(MUST_HAVE_PROU_ERR)
    elif (mod, deg) not in CACHED_FIND_PRIMITIVE_ROOT:
        CACHED_FIND_PRIMITIVE_ROOT[(mod, deg)] = None
        r: int = 2
        while r < mod and not is_prou(root=r, mod=mod, deg=deg):
            r += 1
        if not is_prou(root=r, mod=mod, deg=deg):
            raise RuntimeError(NO_PROU_FOUND_ERR)
        CACHED_FIND_PRIMITIVE_ROOT[(mod, deg)] = r
    return CACHED_FIND_PRIMITIVE_ROOT[(mod, deg)]


def is_root_inverse(root: int, inv_root: int, mod: int) -> bool:
    return (root * inv_root) % mod == 1


def check_ntt_and_intt(val: List[int], mod: int, halfmod: int, logmod: int, deg: int, root: int, inv_root: int,
                       inv_flag: bool, brv_powers: List[int]):
    if not isinstance(val, list) or not isinstance(brv_powers, list):
        raise TypeError(MUST_BE_LIST_ERR)
    elif not all(isinstance(x, int) for x in val + brv_powers + [mod, deg, halfmod, logmod]):
        raise TypeError(MUST_BE_INT_ERR)
    elif not isinstance(inv_flag, bool):
        raise TypeError(MUST_BE_BOOL_ERR)
    elif not is_odd_prime(val=mod):
        raise ValueError(MUST_BE_ODD_PRIME_ERR)
    elif not has_prou(mod=mod, deg=deg):
        raise ValueError(MUST_HAVE_PROU_ERR)
    elif not is_pow_two_geq_two(val=len(val)):
        raise ValueError(MUST_BE_LIST_W_POW_2_LEN_ERR)
    elif deg != len(val):
        raise ValueError(DEGREE_MISMATCH_ERR)
    elif not is_prou(root=root, mod=mod, deg=deg):
        raise ValueError(MUST_BE_CORRECT_ROOT_ERR)
    elif not is_root_inverse(root=root, inv_root=inv_root, mod=mod):
        raise ValueError(MUST_BE_CORRECT_INVERSE_ROOT_ERR)
    elif (not inv_flag and bit_reverse_copy([pow(base=root, exp=i, mod=mod) for i in range(deg)]) != brv_powers) or (inv_flag and bit_reverse_copy([pow(base=inv_root, exp=i, mod=mod) for i in range(deg)]) != brv_powers):
        raise ValueError(MUST_CONSTRUCT_BRV_POWERS_CORRECTLY_ERR)
    check_mod_halfmod_logmod(mod=mod, halfmod=halfmod, logmod=logmod)
    # not isinstance(mod, int) or not all(isinstance(x, int) for x in val) or not isinstance(deg, int) or
    # not all(isinstance(x, int) for x in brv_powers) or not isinstance(halfmod, int) or not isinstance(logmod, int):


def cent(val: int, mod: int, halfmod: int, logmod: int) -> int:
    """
    Centrally reduce a value modulo q in constant time. Output val satisfies -(q//2) <= val <= q//2.
    """
    check_mod_halfmod_logmod(mod=mod, halfmod=halfmod, logmod=logmod)
    y: int = val % mod
    w: int = y - halfmod - 1
    z: int = y - (1 + (w >> logmod)) * mod
    return z


def cooley_tukey_ntt(val: List[int], mod: int, halfmod: int, logmod: int, deg: int, root, brv_powers: List[int]) -> List[int]:
    """
    Input val, a list of n := len(val) integers in usual ordering, a modulus that is a prime such that
    (modulus-1) % (2*n) == 0, a root_order == 2*n, and a list of integers, root_powers, which are powers of a primitive
    root of unity with order root_order stored in bit-reversed order.  Output the NTT of val in bit-reversed order.

    In-place computation, iterative implementation of the Cooley-Tukey butterfly as defined in 'Speeding up the Number
    Theoretic Transform for Faster Ideal Lattice-Based Cryptography' by Longa and Naehrig
    (https://eprint.iacr.org/2016/504.pdf).
    """
    check_ntt_and_intt(val=val, mod=mod, halfmod=halfmod, logmod=logmod, deg=deg, root=root,
                       inv_root=pow(base=root, exp=mod - 2, mod=mod), inv_flag=False, brv_powers=brv_powers)
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
            s: int = brv_powers[m + i]
            for j in range(j_one, j_two + 1):
                u, v = val[j], val[j + t] * s
                val[j] = cent(val=u + v, mod=mod, halfmod=halfmod, logmod=logmod)
                val[j + t] = cent(val=u - v, mod=mod, halfmod=halfmod, logmod=logmod)
        m *= 2
    return val


def gentleman_sande_intt(val: List[int], mod: int, halfmod: int, logmod: int, deg: int, inv_root,
                         brv_powers: List[int]) -> List[int]:
    """
    Input val, a list of n := len(val) integers in usual ordering, a modulus that is a prime such that
    (modulus-1) % (2*n) == 0, a root_order == 2*n, and a list of integers, inv_root_powers, which are powers of the
    inverse of the primitive root of unity with order root_order used to compute the forward transform in bit-reversed
    order.  Output the INTT of val in standard order.

    In-place computation, iterative implementation of the Gentleman-Sande butterfly as in 'Speeding up the Number
    Theoretic Transform for Faster Ideal Lattice-Based Cryptography' by Longa and Naehrig
    (https://eprint.iacr.org/2016/504.pdf).
    """
    check_ntt_and_intt(val=val, mod=mod, halfmod=halfmod, logmod=logmod, deg=deg,
                       root=pow(base=inv_root, exp=mod - 2, mod=mod), inv_root=inv_root, inv_flag=True,
                       brv_powers=brv_powers)
    u: int
    v: int
    n: int = len(val)
    n_inv: int = pow(n, mod - 2, mod)
    t: int = 1
    m: int = n
    while m > 1:
        j_one: int = 0
        h: int = m // 2
        for i in range(h):
            j_two: int = j_one + t - 1
            s: int = brv_powers[h + i]
            for j in range(j_one, j_two + 1):
                u, v = val[j], val[j + t]
                val[j] = cent(val=u + v, mod=mod, halfmod=halfmod, logmod=logmod)
                val[j + t] = cent(val=(u - v) * s, mod=mod, halfmod=halfmod, logmod=logmod)
            j_one += 2 * t
        t *= 2
        m //= 2
    for j in range(n):
        val[j] = cent(val=val[j] * n_inv, mod=mod, halfmod=halfmod, logmod=logmod)
    return val


def check_ntt_poly_mult(f: List[int], g: List[int], deg: int, mod: int, halfmod: int, logmod: int, root_order: int,
                        root: int, inv_root: int, brv_powers: List[int]):
    # Re-use existing check functions for common checks
    check_mod_halfmod_logmod(mod, halfmod, logmod)

    # Check for both f and g using check_ntt_and_intt
    for val in [f, g]:
        check_ntt_and_intt(val=val, mod=mod, halfmod=halfmod, logmod=logmod, deg=deg, root=root,
                           inv_root=inv_root, inv_flag=False, brv_powers=brv_powers)

    # Perform unique checks for ntt_poly_mult function
    if not isinstance(root, int) or not isinstance(inv_root, int):
        raise TypeError(MUST_BE_INT_ERR)
    elif not has_prou(mod=mod, deg=deg):
        raise ValueError(MUST_HAVE_PROU_ERR)
    elif not is_prou(root=root, mod=mod, deg=deg):
        raise ValueError(MUST_BE_CORRECT_ROOT_ERR)
    elif not is_root_inverse(root=root, inv_root=inv_root, mod=mod):
        raise ValueError(MUST_BE_CORRECT_INVERSE_ROOT_ERR)
    elif len(f) != root_order//2:
        raise ValueError(MUST_BE_HALF_FLOORED_ERR)
    elif len(f) != len(g) or len(f) != deg:
        raise ValueError(DEGREE_MISMATCH_ERR)
    pass


def ntt_poly_mult(f: List[int], g: List[int], mod: int, halfmod: int, logmod: int, deg: int, root_order: int, root: int,
                  inv_root: int, brv_powers: List[int], brv_inv_root_powers: List[int]) -> List[int]:
    """
    Multiply two coefficient representations of polynomials by first computing their NTTs
    and then their Hadamard product before inverting the NTT.
    """
    check_ntt_poly_mult(f=f, g=g, deg=deg, mod=mod, halfmod=halfmod, logmod=logmod, root_order=root_order, root=root, inv_root=inv_root, brv_powers=brv_powers)
    cooley_tukey_ntt(val=f, mod=mod, halfmod=halfmod, logmod=logmod, deg=deg, root=root, brv_powers=brv_powers)
    cooley_tukey_ntt(val=g, mod=mod, halfmod=halfmod, logmod=logmod, deg=deg, root=root, brv_powers=brv_powers)
    fg: List[int] = [cent(val=x * y, mod=mod, halfmod=halfmod, logmod=logmod) for x, y in zip(f, g)]
    gentleman_sande_intt(val=fg, mod=mod, halfmod=halfmod, logmod=logmod, deg=deg, inv_root=inv_root, brv_powers=brv_inv_root_powers)
    gentleman_sande_intt(val=f, mod=mod, halfmod=halfmod, logmod=logmod, deg=deg, inv_root=inv_root, brv_powers=brv_inv_root_powers)
    gentleman_sande_intt(val=g, mod=mod, halfmod=halfmod, logmod=logmod, deg=deg, inv_root=inv_root, brv_powers=brv_inv_root_powers)
    return fg
