"""
The ntt module handles the Number Theoretic Transform (NTT) and its inverse in constant time.
"""
from math import ceil
from typing import Dict, Tuple, List, Union
from copy import deepcopy

CACHED_PRIMITIVE_ROOTS: Dict[Tuple[int, int], int] = {}
CACHED_IS_ODD_PRIME: Dict[int, bool] = {}
CACHED_HAS_PRIMITIVE_ROOT_OF_UNITY: Dict[Tuple[int, int], bool] = {}
CACHED_IS_POW_TWO_GEQ_TWO: Dict[int, bool] = {}
CACHED_IS_ROOT_OF_UNITY: Dict[Tuple[int, int, int], bool] = {}
CACHED_IS_PRIMITIVE_ROOT_OF_UNITY: Dict[Tuple[int, int, int], bool] = {}
CACHED_FIND_PRIMITIVE_ROOT: Dict[Tuple[int, int], Union[int, None]] = {}


def is_odd_prime(val: int) -> bool:
    """
    Check if x is an odd prime number.
    :param val: Input number.
    :type val: int
    :return b: Boolean indicating whether x is an odd prime.
    :rtype: bool
    """
    if not isinstance(val, int):
        raise TypeError(f"Value must be an integer, not {type(val)}")
    elif val not in CACHED_IS_ODD_PRIME:
        CACHED_IS_ODD_PRIME[val] = val >= 3 and all(val % i != 0 for i in range(2, ceil(val**0.5) + 1))
    return CACHED_IS_ODD_PRIME[val]


def has_primitive_root_of_unity(modulus: int, root_order: int) -> bool:
    """
    Check if a modulus has a primitive root of unity of a given order.
    :param modulus: Integer modulus
    :type modulus: int
    :param root_order: Order of primitive root of unity
    :type root_order: int
    :return: True if modulus has a primitive root of unity of order root_order, False otherwise
    """
    if not isinstance(modulus, int) or not isinstance(root_order, int):
        raise TypeError(f"Modulus and root order must be integers, not {(type(modulus), type(root_order))}")
    elif (modulus, root_order) not in CACHED_HAS_PRIMITIVE_ROOT_OF_UNITY:
        CACHED_HAS_PRIMITIVE_ROOT_OF_UNITY[(modulus, root_order)] = (
                modulus >= 3
                and root_order >= 2
                and (modulus - 1) % root_order == 0)
    return CACHED_HAS_PRIMITIVE_ROOT_OF_UNITY[(modulus, root_order)]


def is_pow_two_geq_two(val: int) -> bool:
    """
    Check if a number is a power of two.
    :param val: Input integer
    :type val: int
    :return: True if val is a power of two, False otherwise
    :rtype: bool
    """
    if not isinstance(val, int):
        raise TypeError(f"val must be an integer")
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


def check_modulus_halfmod_logmod(modulus: int, halfmod: int, logmod: int):
    if any(not isinstance(x, int) for x in [modulus, halfmod, logmod]):
        raise TypeError("Modulus, halfmod, and logmod must all be integers.")
    elif modulus < 3:
        raise ValueError(f"Modulus must be >=3, had {modulus}")
    elif halfmod < 1:
        raise ValueError(f"Halfmod must be >=1, had {halfmod}")
    elif halfmod != modulus//2:
        raise ValueError(f"Halfmod must be half modulus, but had halfmod={halfmod} and modulus={modulus}")
    elif logmod < 1:
        raise ValueError(f"Logmod must be >=1, had {logmod}")
    elif 2**(logmod-1) >= modulus or modulus > 2**logmod:
        raise ValueError(f"Logmod must bound the modulus, had 2**(logmod-1)={2**(logmod-1)}, 2**logmod={2**logmod},"
                         + f" and modulus={modulus}")
    pass


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
    check_modulus_halfmod_logmod(modulus=modulus, halfmod=halfmod, logmod=logmod)
    y: int = val % modulus
    intermediate_value: int = y - halfmod - 1
    z: int = y - (1 + (intermediate_value >> logmod)) * modulus
    return z


def is_root_of_unity(purported_root: int, modulus: int, root_order: int) -> bool:
    """
    Check if val is a root of unity of order root_order modulo modulus.
    :param purported_root: Input integer
    :type purported_root: int
    :param modulus: Input modulus
    :type modulus: int
    :param root_order: Order of the root of unity
    :type root_order: int
    :return b: Boolean indicating whether val is a root of unity of order root_order modulo modulus
    :rtype: bool
    """
    if not isinstance(purported_root, int) or not isinstance(modulus, int) or not isinstance(root_order, int):
        raise TypeError("Input must be integers")
    if (purported_root, modulus, root_order) not in CACHED_IS_ROOT_OF_UNITY:
        CACHED_IS_ROOT_OF_UNITY[(purported_root, modulus, root_order)] = (
                modulus >= 3
                and root_order >= 2
                and pow(purported_root, root_order, modulus) == 1)
    return CACHED_IS_ROOT_OF_UNITY[(purported_root, modulus, root_order)]


def is_primitive_root(purported_root: int, modulus: int, root_order: int) -> bool:
    """
    Check if val is a primitive root of order root_order modulo modulus.
    :param purported_root: Input integer
    :type purported_root: int
    :param modulus: Input modulus
    :type modulus: int
    :param root_order: Order of the root of unity
    :type root_order: int
    :return b: Boolean indicating whether val is a primitive root of order root_order modulo modulus
    :rtype: bool
    """
    if not isinstance(purported_root, int) or not isinstance(modulus, int) or not isinstance(root_order, int):
        raise TypeError("Input must be integers")
    if (purported_root, modulus, root_order) not in CACHED_IS_PRIMITIVE_ROOT_OF_UNITY:
        is_rou: bool = is_root_of_unity(purported_root=purported_root, modulus=modulus, root_order=root_order)
        is_prim: bool = is_rou and all(pow(purported_root, i, modulus) != 1 for i in range(1, root_order))
        CACHED_IS_PRIMITIVE_ROOT_OF_UNITY[(purported_root, modulus, root_order)] = is_prim
    return CACHED_IS_PRIMITIVE_ROOT_OF_UNITY[(purported_root, modulus, root_order)]


def find_primitive_root(modulus: int, root_order: int) -> int:
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
        raise TypeError("modulus and root_order must be integers")
    elif not has_primitive_root_of_unity(modulus=modulus, root_order=root_order):
        raise ValueError("modulus does not have primitive root of unity of correct order")
    elif (modulus, root_order) not in CACHED_FIND_PRIMITIVE_ROOT:
        CACHED_FIND_PRIMITIVE_ROOT[(modulus, root_order)] = None
        r: int = 2
        while r < modulus and not is_primitive_root(purported_root=r, modulus=modulus, root_order=root_order):
            r += 1
        if not is_primitive_root(purported_root=r, modulus=modulus, root_order=root_order):
            raise RuntimeError(
                f"No primitive root found with modulus={modulus}, root_order={root_order}, despite that"
                + f"has_primitive_root_of_unity returned True."
            )
        CACHED_FIND_PRIMITIVE_ROOT[(modulus, root_order)] = r
    return CACHED_FIND_PRIMITIVE_ROOT[(modulus, root_order)]


def check_ntt_and_intt(val: List[int], modulus: int, root_order: int, brv_powers: List[int], inv_flag: bool,
                       halfmod: int, logmod: int):
    if not isinstance(val, list) or not all(isinstance(x, int) for x in val):
        raise TypeError("input val must be a list of integers")
    elif not isinstance(modulus, int):
        raise TypeError("modulus must be an integer")
    elif not isinstance(root_order, int):
        raise TypeError("root_order must be an integer")
    elif not isinstance(brv_powers, list) or not all(isinstance(x, int) for x in brv_powers):
        raise TypeError("brv_powers must be a list of integers")
    elif not isinstance(inv_flag, bool):
        raise TypeError("inv_flag must be an boolean")
    elif not isinstance(halfmod, int):
        raise TypeError("halfmod must be an integer")
    elif not isinstance(logmod, int):
        raise TypeError("logmod must be a integer")
    elif not is_odd_prime(val=modulus):
        raise ValueError(f"modulus={modulus} must be an odd prime")
    elif not has_primitive_root_of_unity(modulus=modulus, root_order=root_order):
        raise ValueError(
            f"modulus={modulus} does not have a primitive root of order root_order={root_order}"
        )
    elif not is_pow_two_geq_two(val=len(val)):
        raise ValueError(f"len(val)={len(val)} must be a power of 2 greater than 1")
    elif root_order != 2 * len(val):
        raise ValueError(
            f"root_order={root_order} must be twice the degree, {len(val)}"
        )
    elif halfmod != modulus//2:
        raise ValueError(f"halfmod must be half the modulus")
    elif 2**(logmod-1) >= modulus or 2**logmod < modulus:
        raise ValueError(f"logmod does not properly capture bit length of modulus")
    elif (not inv_flag
          and find_primitive_root(modulus=modulus, root_order=root_order) != bit_reverse_copy(brv_powers)[1]):
        fpr = find_primitive_root(modulus=modulus, root_order=root_order)
        raise ValueError(
            f"Need root(={fpr}) to == brv_power = {bit_reverse_copy(brv_powers)[1]} to compute forward NTT"
        )
    elif (pow(
            find_primitive_root(modulus=modulus, root_order=root_order),
            modulus - 2,
            modulus)
          != bit_reverse_copy(brv_powers)[1] and inv_flag):
        raise ValueError(
            f"Need inv_root={pow(find_primitive_root(modulus=modulus, root_order=root_order), modulus - 2, modulus)}"
            + f" == {bit_reverse_copy(brv_powers)[1]} for INTT to be valid"
        )
    pass


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
    check_ntt_and_intt(val=val, modulus=modulus, root_order=root_order, brv_powers=bit_rev_root_powers, inv_flag=False,
                       halfmod=halfmod, logmod=logmod)
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
                val[j] = cent(val=u + v, modulus=modulus, halfmod=halfmod, logmod=logmod)
                val[j + t] = cent(val=u - v, modulus=modulus, halfmod=halfmod, logmod=logmod)
        m *= 2
    return val


def gentleman_sande_intt(val: List[int], modulus: int, root_order: int, bit_rev_inv_root_powers: List[int],
                         halfmod: int, logmod: int) -> List[int]:
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
    check_ntt_and_intt(val=val, modulus=modulus, root_order=root_order, brv_powers=bit_rev_inv_root_powers,
                       inv_flag=True, halfmod=halfmod, logmod=logmod)
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
        m = m // 2
    for j in range(n):
        val[j] = cent(val=val[j] * n_inv, modulus=modulus, halfmod=halfmod, logmod=logmod)
    return val


def check_ntt_poly_mult(f: List[int], g: List[int], degree: int, modulus: int, halfmod: int, logmod: int,
                        root_order: int, root: int, inv_root: int, brv_root_powers: List[int]):
    # Re-use existing check functions for common checks
    check_modulus_halfmod_logmod(modulus, halfmod, logmod)

    # Check for both f and g using check_ntt_and_intt
    for val in [f, g]:
        check_ntt_and_intt(val=val, modulus=modulus, root_order=root_order, brv_powers=brv_root_powers, inv_flag=False,
                           halfmod=halfmod, logmod=logmod)

    # Perform unique checks for ntt_poly_mult function
    if not isinstance(root, int) or not isinstance(inv_root, int):
        raise TypeError("Parameters root and inv_root must be integers.")
    elif not is_primitive_root(purported_root=root, modulus=modulus, root_order=root_order):
        raise ValueError("Root must be a primitive root of unity.")
    elif (root * inv_root) % modulus != 1:
        raise ValueError("Inv_root must be the inverse of the root of unity.")
    elif len(f) != len(g) or len(f) != root_order // 2 or len(f) != degree:
        raise ValueError("f and g must have the same length and be of length root_order // 2.")
    pass


def ntt_poly_mult(f: List[int], g: List[int], modulus: int, halfmod: int, logmod: int, degree: int, root_order: int,
                  root: int, inv_root: int, brv_root_powers: List[int]) -> List[int]:
    """
    Multiply two coefficient representations of polynomials by first computing their NTTs
    and then their Hadamard product before inverting the NTT.
    """
    check_ntt_poly_mult(f=f, g=g, modulus=modulus, halfmod=halfmod, logmod=logmod, degree=degree, root_order=root_order,
                        root=root, inv_root=inv_root, brv_root_powers=brv_root_powers)
    root_powers: List[int] = [pow(root, i, modulus) for i in range(degree)]
    brv_root_powers: List[int] = bit_reverse_copy(val=root_powers)
    inv_root_powers: List[int] = [pow(inv_root, i, modulus) for i in range(degree)]
    brv_inv_root_powers = bit_reverse_copy(val=inv_root_powers)
    fhat = cooley_tukey_ntt(val=deepcopy(f), modulus=modulus, root_order=root_order, bit_rev_root_powers=brv_root_powers,
                     halfmod=halfmod, logmod=logmod)
    ghat = cooley_tukey_ntt(val=deepcopy(g), modulus=modulus, root_order=root_order, bit_rev_root_powers=brv_root_powers,
                     halfmod=halfmod, logmod=logmod)
    fghat: List[int] = [
        cent(val=x * y, modulus=modulus, halfmod=halfmod, logmod=logmod)
        for x, y in zip(fhat, ghat)
    ]
    return gentleman_sande_intt(val=fghat, modulus=modulus, root_order=root_order, bit_rev_inv_root_powers=brv_inv_root_powers,
                         halfmod=halfmod, logmod=logmod)


def derived_params(val: List[int], modulus: int, inv_flag: bool) -> Tuple[int, int, int, int, int, List[int]]:
    """ Derive root order from list length, derive root or inv_root from modulus and root order and inv_flag,
    compute root powers, then bit reverse."""
    if not isinstance(val, list) or not all(isinstance(x, int) for x in val) or not isinstance(modulus, int):
        raise TypeError(f"Input must be lists of integers and integer.")
    elif not is_odd_prime(val=modulus):
        raise ValueError(f"Modulus must be an odd prime")
    elif not has_primitive_root_of_unity(modulus=modulus, root_order=2*len(val)):
        raise ValueError(f"Modulus does not have primitive root of unity of correct order")
    degree: int = len(val)
    root_order: int = 2*degree
    halfmod: int = modulus//2
    logmod: int = modulus.bit_length()
    root: int = find_primitive_root(modulus=modulus, root_order=root_order)
    if inv_flag:
        root = pow(base=root, exp=modulus - 2, mod=modulus)
    powers: List[int] = [pow(base=root, exp=i, mod=modulus) for i in range(degree)]
    brv_powers: List[int] = bit_reverse_copy(powers)
    return halfmod, logmod, degree, root_order, root, brv_powers,
