"""
The ntt module handles the Number Theoretic Transform (NTT) and its inverse in constant time.
"""
from math import ceil
from typing import Dict, Tuple, List, Union
from copy import deepcopy

# Caching dictionaries

_CACHED_PRIMITIVE_ROOTS: Dict[Tuple[int, int], int] = {}
_CACHED_IS_ODD_PRIME: Dict[int, bool] = {}
_CACHED_HAS_PRIMITIVE_ROOT_OF_UNITY: Dict[Tuple[int, int], bool] = {}
_CACHED_IS_POW_TWO_GEQ_TWO: Dict[int, bool] = {}
_CACHED_IS_ROOT_OF_UNITY: Dict[Tuple[int, int, int], bool] = {}
_CACHED_IS_PRIMITIVE_ROOT_OF_UNITY: Dict[Tuple[int, int, int], bool] = {}
_CACHED_FIND_PRIMITIVE_ROOT: Dict[Tuple[int, int], Union[int, None]] = {}

# Error messages

_ERR_MUST_BE_INT: str = 'val must be an int'
_ERR_MUST_BE_LIST: str = 'val must be a list'
_ERR_MUST_BE_BOOL: str = 'val must be a bool'
_ERR_MUST_BE_POW_TWO: str = 'val must be >= 2 and power of 2'
_ERR_MUST_BE_GEQ_THREE: str = 'val must be >= 3'
_ERR_MUST_BE_STRICTLY_POSITIVE: str = 'val must be > 0'
_ERR_MUST_BE_FLOOR_HALF: str = 'val must be floor-divided by 2'
_ERR_MUST_BE_BIT_LENGTH: str = 'val must be correct bit-length'
_ERR_MUST_HAVE_PRIM_ROU: str = 'val must have primitive root of unity of correct order'
_ERR_MUST_BE_ODD_PRIME: str = 'val must be odd prime'
_ERR_MUST_BE_CORRECT_ROOT: str = 'val must be correct root'
_ERR_MUST_BE_SAME_LENGTH: str = 'vals must be same length'


def _is_odd_prime(val: int) -> bool:
    """
    Check if a number is an odd prime.
    :param val: Input integer
    :type val: int
    :return: True if val is an odd prime, False otherwise
    """
    if not isinstance(val, int):
        raise TypeError(_ERR_MUST_BE_INT)
    elif val not in _CACHED_IS_ODD_PRIME:
        _CACHED_IS_ODD_PRIME[val] = val >= 3 and all(val % i != 0 for i in range(2, ceil(val ** 0.5) + 1))
    return _CACHED_IS_ODD_PRIME[val]


def _has_primitive_root_of_unity(modulus: int, root_order: int) -> bool:
    """
    Check if a modulus has a primitive root of unity of a given order.
    :param modulus: Integer modulus
    :type modulus: int
    :param root_order: Order of primitive root of unity
    :type root_order: int
    :return: True if modulus has a primitive root of unity of order root_order, False otherwise
    """
    if not isinstance(modulus, int) or not isinstance(root_order, int):
        raise TypeError(_ERR_MUST_BE_INT)
    elif (modulus, root_order) not in _CACHED_HAS_PRIMITIVE_ROOT_OF_UNITY:
        _CACHED_HAS_PRIMITIVE_ROOT_OF_UNITY[(modulus, root_order)] = (
                modulus >= 3
                and root_order >= 2
                and (modulus - 1) % root_order == 0)
    return _CACHED_HAS_PRIMITIVE_ROOT_OF_UNITY[(modulus, root_order)]


def _is_pow_two_geq_two(val: int) -> bool:
    """
    Check if a number is a power of two.
    :param val: Input integer
    :type val: int
    :return: True if val is a power of two, False otherwise
    :rtype: bool
    """
    if not isinstance(val, int):
        raise TypeError(_ERR_MUST_BE_INT)
    elif val not in _CACHED_IS_POW_TWO_GEQ_TWO:
        _CACHED_IS_POW_TWO_GEQ_TWO[val] = val >= 2 and val & (val - 1) == 0
    return _CACHED_IS_POW_TWO_GEQ_TWO[val]


def _bit_reverse_copy(val: List[int]):
    """
    Permute indices by bit-reversal.
    :param val: Input list of anything
    :type val: list
    :return lav: Output list matches input list with indices bit-reversed
    :rtype: list
    """
    if not isinstance(val, list):
        raise TypeError(_ERR_MUST_BE_LIST)
    elif not all(isinstance(x, int) for x in val):
        raise TypeError(_ERR_MUST_BE_INT)
    elif not _is_pow_two_geq_two(val=len(val)):
        raise ValueError(_ERR_MUST_BE_POW_TWO)
    n: int = len(val)
    k: int = n.bit_length() - 1
    bit_reversed_indices: List[int] = [int(bin(i)[2:].zfill(k)[::-1], 2) for i in range(n)]
    lav = [val[i] for i in bit_reversed_indices]  # bit reverse
    return lav


def _check_modulus_halfmod_logmod(modulus: int, halfmod: int, logmod: int):
    if any(not isinstance(x, int) for x in [modulus, halfmod, logmod]):
        raise TypeError(_ERR_MUST_BE_INT)
    elif halfmod < 1 or logmod < 1:
        raise ValueError(_ERR_MUST_BE_STRICTLY_POSITIVE)
    elif modulus < 3:
        raise ValueError(_ERR_MUST_BE_GEQ_THREE)
    elif halfmod != modulus//2:
        raise ValueError(_ERR_MUST_BE_FLOOR_HALF)
    elif modulus.bit_length() != logmod:
        raise ValueError(_ERR_MUST_BE_BIT_LENGTH)
    pass


def _cent(val: int, modulus: int, halfmod: int, logmod: int) -> int:
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
    _check_modulus_halfmod_logmod(modulus=modulus, halfmod=halfmod, logmod=logmod)
    y: int = val % modulus
    intermediate_value: int = y - halfmod - 1
    z: int = y - (1 + (intermediate_value >> logmod)) * modulus
    return z


def _is_root_of_unity(val: int, modulus: int, root_order: int) -> bool:
    """
    Check if val is a root of unity of order root_order modulo modulus.
    :param val: Input integer
    :type val: int
    :param modulus: Input modulus
    :type modulus: int
    :param root_order: Order of the root of unity
    :type root_order: int
    :return b: Boolean indicating whether val is a root of unity of order root_order modulo modulus
    :rtype: bool
    """
    if not isinstance(val, int) or not isinstance(modulus, int) or not isinstance(root_order, int):
        raise TypeError(_ERR_MUST_BE_INT)
    if (val, modulus, root_order) not in _CACHED_IS_ROOT_OF_UNITY:
        _CACHED_IS_ROOT_OF_UNITY[(val, modulus, root_order)] = (
                modulus >= 3
                and root_order >= 2
                and pow(val, root_order, modulus) == 1)
    return _CACHED_IS_ROOT_OF_UNITY[(val, modulus, root_order)]


def _is_primitive_root(val: int, modulus: int, root_order: int) -> bool:
    """
    Check if val is a primitive root of order root_order modulo modulus.
    :param val: Input integer
    :type val: int
    :param modulus: Input modulus
    :type modulus: int
    :param root_order: Order of the root of unity
    :type root_order: int
    :return b: Boolean indicating whether val is a primitive root of order root_order modulo modulus
    :rtype: bool
    """
    if not isinstance(val, int) or not isinstance(modulus, int) or not isinstance(root_order, int):
        raise TypeError(_ERR_MUST_BE_INT)
    if (val, modulus, root_order) not in _CACHED_IS_PRIMITIVE_ROOT_OF_UNITY:
        is_rou: bool = _is_root_of_unity(val=val, modulus=modulus, root_order=root_order)
        is_prim: bool = is_rou and all(pow(val, i, modulus) != 1 for i in range(1, root_order))
        _CACHED_IS_PRIMITIVE_ROOT_OF_UNITY[(val, modulus, root_order)] = is_prim
    return _CACHED_IS_PRIMITIVE_ROOT_OF_UNITY[(val, modulus, root_order)]


def _find_primitive_root(modulus: int, root_order: int) -> int:
    """
    Find a primitive root of order root_order modulo modulus. Naive loop that first checks 2, then 3, then 4...
    :param modulus: Input modulus
    :type modulus: int
    :param root_order: Order of the root of unity
    :type root_order: int
    :return r: root of unity r such that r**root_order % q == 1 and r**i % q != 1 for all i in range(1, root_order)
    :rtype: int
    """
    if not isinstance(modulus, int) or not isinstance(root_order, int):
        raise TypeError(_ERR_MUST_BE_INT)
    elif not _has_primitive_root_of_unity(modulus=modulus, root_order=root_order):
        raise ValueError(_ERR_MUST_HAVE_PRIM_ROU)
    elif (modulus, root_order) not in _CACHED_FIND_PRIMITIVE_ROOT:
        _CACHED_FIND_PRIMITIVE_ROOT[(modulus, root_order)] = None
        r: int = 2
        while r < modulus and not _is_primitive_root(val=r, modulus=modulus, root_order=root_order):
            r += 1
        if not _is_primitive_root(val=r, modulus=modulus, root_order=root_order):
            raise RuntimeError(_ERR_MUST_HAVE_PRIM_ROU)
        _CACHED_FIND_PRIMITIVE_ROOT[(modulus, root_order)] = r
    return _CACHED_FIND_PRIMITIVE_ROOT[(modulus, root_order)]


def _check_ntt_and_intt(val: List[int], modulus: int, halfmod: int, logmod: int, root_order: int, brv_powers: List[int], inv_flag: bool):
    if not isinstance(val, list) or not isinstance(brv_powers, list):
        raise TypeError(_ERR_MUST_BE_LIST)
    elif not all(isinstance(x, int) for x in val) or not isinstance(modulus, int) or not isinstance(root_order, int) or not all(isinstance(x, int) for x in brv_powers) or not isinstance(halfmod, int) or not isinstance(logmod, int):
        raise TypeError(_ERR_MUST_BE_INT)
    elif not isinstance(inv_flag, bool):
        raise TypeError(_ERR_MUST_BE_BOOL)
    elif not _is_odd_prime(val=modulus):
        raise ValueError(_ERR_MUST_BE_ODD_PRIME)
    elif not _has_primitive_root_of_unity(modulus=modulus, root_order=root_order):
        raise ValueError(_ERR_MUST_HAVE_PRIM_ROU)
    elif not _is_pow_two_geq_two(val=len(val)):
        raise ValueError(_ERR_MUST_BE_POW_TWO)
    elif root_order//2 != len(val) or halfmod != modulus//2:
        raise ValueError(_ERR_MUST_BE_FLOOR_HALF)
    elif modulus.bit_length() != logmod:
        raise ValueError(_ERR_MUST_BE_BIT_LENGTH)
    elif not inv_flag and _find_primitive_root(modulus=modulus, root_order=root_order) != _bit_reverse_copy(brv_powers)[1]:
        fpr = _find_primitive_root(modulus=modulus, root_order=root_order)
        raise ValueError(_ERR_MUST_BE_CORRECT_ROOT)
    elif inv_flag and pow(base=_find_primitive_root(modulus=modulus, root_order=root_order), exp=modulus - 2, mod=modulus) != _bit_reverse_copy(brv_powers)[1]:
        raise ValueError(_ERR_MUST_BE_CORRECT_ROOT)
    pass


def _cooley_tukey_ntt(val: List[int], modulus: int, halfmod: int, logmod: int, root_order: int,
                      bit_rev_root_powers: List[int]) -> List[int]:
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
    _check_ntt_and_intt(val=val, modulus=modulus, halfmod=halfmod, logmod=logmod, root_order=root_order,
                        brv_powers=bit_rev_root_powers, inv_flag=False)
    u: int
    v: int
    n: int = len(val)
    t: int = n
    m: int = 1
    while m < n:
        t //= 2
        for i in range(m):
            j_one: int = 2 * i * t
            j_two: int = j_one + t - 1
            s: int = bit_rev_root_powers[m + i]
            for j in range(j_one, j_two + 1):
                u, v = val[j], val[j + t] * s
                val[j] = _cent(val=u + v, modulus=modulus, halfmod=halfmod, logmod=logmod)
                val[j + t] = _cent(val=u - v, modulus=modulus, halfmod=halfmod, logmod=logmod)
        m *= 2
    return val


def _gentleman_sande_intt(val: List[int], modulus: int, halfmod: int, logmod: int, root_order: int,
                          bit_rev_inv_root_powers: List[int]) -> List[int]:
    """
    Input val, a list of n := len(val) integers in bit-reversed ordering, a modulus that is a prime such that
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
    _check_ntt_and_intt(val=val, modulus=modulus, halfmod=halfmod, logmod=logmod, root_order=root_order,
                        brv_powers=bit_rev_inv_root_powers, inv_flag=True)
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
                val[j] = _cent(val=u + v, modulus=modulus, halfmod=halfmod, logmod=logmod)
                val[j + t] = _cent(val=(u - v) * s, modulus=modulus, halfmod=halfmod, logmod=logmod)
            j_one += 2 * t
        t *= 2
        m = m // 2
    for j in range(n):
        val[j] = _cent(val=val[j] * n_inv, modulus=modulus, halfmod=halfmod, logmod=logmod)
    return val


def _check_ntt_poly_mult(f: List[int], g: List[int], modulus: int, halfmod: int, logmod: int, degree: int,
                         root_order: int, root: int, inv_root: int, brv_root_powers: List[int]):
    # Re-use existing check functions for common checks
    _check_modulus_halfmod_logmod(modulus, halfmod, logmod)

    # Check for both f and g using check_ntt_and_intt
    for val in [f, g]:
        _check_ntt_and_intt(val=val, modulus=modulus, halfmod=halfmod, logmod=logmod, root_order=root_order,
                            brv_powers=brv_root_powers, inv_flag=False)

    # Perform unique checks for ntt_poly_mult function
    if not isinstance(root, int) or not isinstance(inv_root, int):
        raise TypeError(_ERR_MUST_BE_INT)
    elif not _is_primitive_root(val=root, modulus=modulus, root_order=root_order) or (root * inv_root) % modulus != 1:
        raise ValueError(_ERR_MUST_BE_CORRECT_ROOT)
    elif len(f) != len(g):
        raise ValueError(_ERR_MUST_BE_SAME_LENGTH)
    elif len(f) != root_order // 2 or degree != root_order//2:
        raise ValueError(_ERR_MUST_BE_FLOOR_HALF)
    pass


def _ntt_poly_mult(f: List[int], g: List[int], modulus: int, halfmod: int, logmod: int, degree: int, root_order: int,
                   root: int, inv_root: int, brv_root_powers: List[int]) -> List[int]:
    """
    Multiply two coefficient representations of polynomials by first computing their NTTs
    and then their Hadamard product before inverting the NTT.
    """
    _check_ntt_poly_mult(f=f, g=g, modulus=modulus, halfmod=halfmod, logmod=logmod, degree=degree, root_order=root_order,
                         root=root, inv_root=inv_root, brv_root_powers=brv_root_powers)
    root_powers: List[int] = [pow(root, i, modulus) for i in range(degree)]
    brv_root_powers: List[int] = _bit_reverse_copy(val=root_powers)
    inv_root_powers: List[int] = [pow(inv_root, i, modulus) for i in range(degree)]
    brv_inv_root_powers = _bit_reverse_copy(val=inv_root_powers)
    fhat = _cooley_tukey_ntt(val=deepcopy(f), modulus=modulus, halfmod=halfmod, logmod=logmod, root_order=root_order,
                             bit_rev_root_powers=brv_root_powers)
    ghat = _cooley_tukey_ntt(val=deepcopy(g), modulus=modulus, halfmod=halfmod, logmod=logmod, root_order=root_order,
                             bit_rev_root_powers=brv_root_powers)
    fghat: List[int] = [
        _cent(val=x * y, modulus=modulus, halfmod=halfmod, logmod=logmod)
        for x, y in zip(fhat, ghat)
    ]
    return _gentleman_sande_intt(val=fghat, modulus=modulus, halfmod=halfmod, logmod=logmod, root_order=root_order,
                                 bit_rev_inv_root_powers=brv_inv_root_powers)


def _derived_params(modulus: int, degree: int, inv_flag: bool) -> Tuple[int, int, int, int, int, int, List[int], bool]:
    """ Derive root order from list length, derive root or inv_root from modulus and root order and inv_flag,
    compute root powers, then bit reverse."""
    if not isinstance(degree, int) or not isinstance(modulus, int):
        raise TypeError(_ERR_MUST_BE_INT)
    elif not isinstance(inv_flag, bool):
        raise TypeError(_ERR_MUST_BE_BOOL)
    elif not _is_odd_prime(val=modulus):
        raise ValueError(_ERR_MUST_BE_ODD_PRIME)
    elif not _has_primitive_root_of_unity(modulus=modulus, root_order=2 * degree):
        raise ValueError(_ERR_MUST_HAVE_PRIM_ROU)
    root_order: int = 2*degree
    halfmod: int = modulus//2
    logmod: int = modulus.bit_length()
    root_or_inv_root: int = _find_primitive_root(modulus=modulus, root_order=root_order)
    if inv_flag:
        root_or_inv_root = pow(base=root_or_inv_root, exp=modulus - 2, mod=modulus)
    powers: List[int] = [pow(base=root_or_inv_root, exp=i, mod=modulus) for i in range(degree)]
    brv_powers: List[int] = _bit_reverse_copy(powers)
    return modulus, halfmod, logmod, degree, root_order, root_or_inv_root, brv_powers, inv_flag
