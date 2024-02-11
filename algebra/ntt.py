"""
algebra/ntt.py

Handle the Number Theoretic Transform (NTT) and its inverse in constant time.

Author: Brandon G Goodell
"""
from math import ceil as _ceil
from copy import deepcopy as _deepcopy
from algebra.errors import (_MustBeIntError, _MustBeListError, _MustBePosPowerOfTwoError, _NoPROUFoundError,
                            _MustBeOddPrimeError, _InvalidNTTInputError, _MustHavePROUError, _MustBeBoolError)

# Caching dictionaries
_ODD_PRIME_CACHE: dict[int, bool] = {}
_HAS_PROU_CACHE: dict[tuple[int, int], bool] = {}
_POW_TWO_CACHE: dict[int, bool] = {}
_ROU_CACHE: dict[tuple[int, int, int], bool] = {}
_PROU_CACHE: dict[tuple[int, int, int], bool] = {}
_ROOT_INV_CACHE: dict[tuple[int, int, int], bool] = {}
_FIND_PROU_CACHE: dict[tuple[int, int], None | int] = {}
_REV_IDX_CACHE = {}
_MOD_HALFMOD_LOGMOD_CACHE: dict[int, tuple[int, int]] = {}
_BRV_ROOTS_AND_INV_ROOTS_CACHE: dict[tuple[int, int], tuple[list[int], list[int]]] = {}


# Boolean Functions
def _is_odd_prime(val: int) -> bool:
    """
    Private, cached function for checking if x is an odd prime number.

    :param val: An integer to check for primality.
    :type val: int
    :return: True if val is an odd prime, False otherwise.
    :rtype: bool
    """
    if val not in _ODD_PRIME_CACHE:
        _ODD_PRIME_CACHE[val] = isinstance(val, int)
        _ODD_PRIME_CACHE[val] = _ODD_PRIME_CACHE[val] and val >= 3
        _ODD_PRIME_CACHE[val] = _ODD_PRIME_CACHE[val] and all(val % i != 0 for i in range(2, _ceil(val ** 0.5) + 1))
    return _ODD_PRIME_CACHE[val]


def _is_pos_pow_two(val: int) -> bool:
    """
    Private, cached function for checking if a number is a power of 2.

    :param val: An integer to check for being a power of 2.
    :type val: int
    :return: True if val is a power of 2, False otherwise.
    :rtype: bool
    """
    if val not in _POW_TWO_CACHE:
        _POW_TWO_CACHE[val] = isinstance(val, int)
        _POW_TWO_CACHE[val] = _POW_TWO_CACHE[val] and val >= 1
        _POW_TWO_CACHE[val] = _POW_TWO_CACHE[val] and val & (val - 1) == 0
    return _POW_TWO_CACHE[val]


def _is_ntt_friendly(mod: int, deg: int) -> bool:
    """
    Private, cached function for checking if a modulus has a primitive root of unity of order 2*deg.

    :param mod: A modulus to check for having a primitive root of unity of order 2*deg.
    :type mod: int
    :param deg: The degree of the polynomial.
    :type deg: int
    :return: True if mod has a primitive root of unity of order 2*deg, False otherwise.
    :rtype: bool
    """
    if (mod, deg) not in _HAS_PROU_CACHE:
        _HAS_PROU_CACHE[(mod, deg)] = isinstance(mod, int)
        _HAS_PROU_CACHE[(mod, deg)] = _HAS_PROU_CACHE[(mod, deg)] and isinstance(deg, int)
        _HAS_PROU_CACHE[(mod, deg)] = _HAS_PROU_CACHE[(mod, deg)] and _is_odd_prime(val=mod)
        _HAS_PROU_CACHE[(mod, deg)] = _HAS_PROU_CACHE[(mod, deg)] and _is_pos_pow_two(val=deg)
        _HAS_PROU_CACHE[(mod, deg)] = _HAS_PROU_CACHE[(mod, deg)] and (mod - 1) % (2 * deg) == 0
    return _HAS_PROU_CACHE[(mod, deg)]


def _is_rou(root: int, mod: int, deg: int) -> bool:
    """
    Private, cached function for determining if an integer is a root of unity of order 2*deg modulo mod.

    :param root: A root of unity to check for having order 2*deg modulo mod.
    :type root: int
    :param mod: A modulus to check for having a primitive root of unity of order 2*deg.
    :type mod: int
    :param deg: The degree of the polynomial.
    :type deg: int
    :return: True if root is a root of unity of order 2*deg modulo mod, False otherwise.
    :rtype: bool
    """
    if (root, mod, deg) not in _ROU_CACHE:
        _ROU_CACHE[(root, mod, deg)] = _is_ntt_friendly(mod=mod, deg=deg)
        _ROU_CACHE[(root, mod, deg)] = _ROU_CACHE[(root, mod, deg)] and isinstance(root, int)
        _ROU_CACHE[(root, mod, deg)] = _ROU_CACHE[(root, mod, deg)] and pow(base=root, exp=2 * deg, mod=mod) == 1
    return _ROU_CACHE[(root, mod, deg)]


def _is_prou(root: int, mod: int, deg: int) -> bool:
    """
    Private, cached function for checking if val is a primitive root of order 2*deg modulo modulus.

    :param root: A root of unity to check for having order 2*deg modulo mod.
    :type root: int
    :param mod: A modulus to check for having a primitive root of unity of order 2*deg.
    :type mod: int
    :param deg: The degree of the polynomial.
    :type deg: int
    :return: True if root is a primitive root of order 2*deg modulo modulus, False otherwise.
    :rtype: bool
    """
    if (root, mod, deg) not in _PROU_CACHE:
        _PROU_CACHE[(root, mod, deg)] = all(pow(base=root, exp=i, mod=mod) != 1 for i in range(1, 2*deg))
        _PROU_CACHE[(root, mod, deg)] = _PROU_CACHE[(root, mod, deg)] and _is_rou(root=root, mod=mod, deg=deg)
    return _PROU_CACHE[(root, mod, deg)]


def _is_root_inverse(root: int, inv_root: int, mod: int) -> bool:
    """
    Private, cached function for checking if root and inv_root are multiplicative inverses modulo mod.

    :param root: A root of unity to check for having order 2*deg modulo mod.
    :type root: int
    :param inv_root: The inverse of root modulo mod.
    :type inv_root: int
    :param mod: A modulus to check for having a primitive root of unity of order 2*deg.
    :type mod: int
    :return: True if root and inv_root are multiplicative inverses modulo mod, False otherwise.
    :rtype: bool
    """
    if (root, inv_root, mod) not in _ROOT_INV_CACHE:
        _ROOT_INV_CACHE[(root, inv_root, mod)] = isinstance(root, int) and isinstance(inv_root, int)
        _ROOT_INV_CACHE[(root, inv_root, mod)] = _ROOT_INV_CACHE[(root, inv_root, mod)] and _is_odd_prime(val=mod)
        _ROOT_INV_CACHE[(root, inv_root, mod)] = _ROOT_INV_CACHE[(root, inv_root, mod)] and (root * inv_root) % mod == 1
    return _ROOT_INV_CACHE[(root, inv_root, mod)]


def _is_ntt_valid(val: list[int], mod: int, inv_flag: bool) -> bool:
    """
    Private function for checking if input to NTT is valid.

    :param val: A list of integers.
    :type val: list
    :param mod: A modulus to check for having a primitive root of unity of order 2*deg.
    :type mod: int
    :param inv_flag: A flag indicating forward or inverse.
    :type inv_flag: bool
    :return: True if the input to NTT is valid, False otherwise.
    :rtype: bool
    """
    val_is_list: bool = isinstance(val, list)
    deg: int = len(val)
    deg_mod_is_ntt_friendly: bool = _is_ntt_friendly(mod=mod, deg=deg)
    inv_flag_is_bool: bool = isinstance(inv_flag, bool)
    entries_in_val_are_ints: bool = all(isinstance(i, int) for i in val)
    brv_powers: list[int]
    brv_powers, brv_inv_powers = _brv_root_and_inv_root_powers(deg=deg, mod=mod)
    powers: list[int] = _brv_copy(val=brv_powers)
    inv_powers: list[int] = _brv_copy(val=brv_inv_powers)
    root: int = powers[1]
    inv_root: int = pow(base=root, exp=mod-2, mod=mod)
    has_correct_powers: bool = all(
        (x-y) % mod == 0 for x, y in zip(powers, [pow(base=root, exp=i, mod=mod) for i in range(deg)]))
    has_correct_inv_powers: bool = all(
        (x-y) % mod == 0 for x, y in zip(inv_powers, [pow(base=inv_root, exp=i, mod=mod) for i in range(deg)]))
    root_is_prou: bool = _is_prou(root=root, mod=mod, deg=deg)
    has_correct_root_inverse: bool = _is_root_inverse(root=root, inv_root=inv_root, mod=mod)
    return (val_is_list and deg_mod_is_ntt_friendly and inv_flag_is_bool and entries_in_val_are_ints and
            has_correct_powers and has_correct_inv_powers and root_is_prou and has_correct_root_inverse)


def _is_ntt_poly_mult_valid(f: list[int], g: list[int], mod: int, brv_powers: list[int]):
    """
    Private function for checking if inputs to NTT polynomial multiplication are valid.

    :param f: Coefficient representation of a polynomial (a list of integers).
    :type f: list[int]
    :param g: Coefficient representation of a polynomial (a list of integers).
    :type g: list[int]
    :param mod: A modulus to reduce val by.
    :type mod: int
    :param brv_powers: A list of powers of the primitive root of unity in bit-reversed order.
    :type brv_powers: list[int]
    :return: True if the inputs to NTT polynomial multiplication are valid, False otherwise.
    :rtype: bool
    """
    deg: int = len(f)
    f_is_list = isinstance(f, list)
    g_is_list = isinstance(g, list)
    mod_is_int = isinstance(mod, int)
    brv_powers_is_list = isinstance(brv_powers, list)
    f_and_g_are_ints = all(isinstance(i, int) for i in f+g)
    f_and_g_have_len_deg = len(f) == len(g) == deg
    have_prou = _is_ntt_friendly(mod=mod, deg=deg)
    powers: list[int] = _brv_copy(val=brv_powers)
    root: int = powers[1]
    inv_root: int = pow(base=root, exp=mod - 2, mod=mod)
    have_correct_root = _is_prou(root=root, mod=mod, deg=deg)
    have_correct_inverse_root = _is_root_inverse(root=root, inv_root=inv_root, mod=mod)
    ntt_is_valid_for_f = _is_ntt_valid(val=f, mod=mod, inv_flag=False)
    ntt_is_valid_for_g = _is_ntt_valid(val=g, mod=mod, inv_flag=False)
    return all([f_is_list, g_is_list, mod_is_int, brv_powers_is_list, f_and_g_are_ints, f_and_g_have_len_deg, have_prou,
                have_correct_root, have_correct_inverse_root, ntt_is_valid_for_f, ntt_is_valid_for_g])


def _find_prou(mod: int, deg: int) -> int:
    """
    Public, cached function for finding a primitive root of order 2*deg modulo mod where mod and deg are ntt-friendly.
    Computes with a naive loop that first checks 2, then 3, then 4...

    :param mod: A modulus to check for having a primitive root of unity of order 2*deg.
    :type mod: int
    :param deg: The degree of the polynomial.
    :type deg: int
    :return: A primitive root of order 2*deg modulo mod.
    :rtype: int
    :raises _MustBeIntError: If mod or deg are not integers
    :raises _MustHavePROUError: If mod and deg are not ntt-friendly
    :raises _NoPROUFoundError: If no primitive root of order 2*deg modulo mod is found.
    """
    if not isinstance(mod, int) or not isinstance(deg, int):
        raise _MustBeIntError
    elif not _is_ntt_friendly(mod=mod, deg=deg):
        raise _MustHavePROUError
    elif (mod, deg) not in _FIND_PROU_CACHE:
        _FIND_PROU_CACHE[(mod, deg)] = None
        r: int = 2
        while r < mod and not _is_prou(root=r, mod=mod, deg=deg):
            r += 1
        if not _is_prou(root=r, mod=mod, deg=deg):
            raise _NoPROUFoundError
        _FIND_PROU_CACHE[(mod, deg)] = r
    return _FIND_PROU_CACHE[(mod, deg)]


def _brv_indices(val: int) -> list[int]:
    """
    Private, cached function for computing bit-reversed indices of a list of length val.

    :param val: The length of the list.
    :type val: int
    :return: A list of bit-reversed indices.
    :rtype: list[int]
    :raises _MustBeIntError: If val is not an integer.
    :raises _MustBePosPowerOfTwoError: If val is not a positive power of two.
    """
    if not isinstance(val, int):
        raise _MustBeIntError
    elif not _is_pos_pow_two(val=val):
        raise _MustBePosPowerOfTwoError
    elif val not in _REV_IDX_CACHE:
        k: int = val.bit_length() - 1
        _REV_IDX_CACHE[val] = [int(bin(i)[2:].zfill(k)[::-1], 2) for i in range(val)]
    return _REV_IDX_CACHE[val]


def _brv_copy(val: list) -> list:
    """
    Public function for computing a copy of 'val', permuting indices by bit-reversal. Do not cache.

    :param val: A list of integers.
    :type val: list
    :return: A list of integers permuted by bit-reversal.
    :rtype: list
    :raises _MustBeListError: If val is not a list.
    """
    if not isinstance(val, list):
        raise _MustBeListError
    brvd_indices: list[int] = _brv_indices(val=len(val))
    return [val[i] for i in brvd_indices]


def _brv_root_and_inv_root_powers(deg: int, mod: int) -> tuple[list[int], list[int]]:
    """
    Private, cached function for computing the primitive root of unity of order 2*deg modulo mod and its inverse, and
    returning their powers.

    :param deg: The degree of the polynomial.
    :type deg: int
    :param mod: A modulus to check for having a primitive root of unity of order 2*deg.
    :type mod: int
    :return: A tuple containing the powers of the primitive root of unity and its inverse.
    :rtype: tuple[list[int], list[int]]
    :raises _MustBeIntError: If deg or mod are not integers.
    :raises _MustHavePROUError: If deg and mod are not ntt-friendly.
    """
    if not isinstance(deg, int) or not isinstance(mod, int):
        raise _MustBeIntError
    elif not _is_ntt_friendly(mod=mod, deg=deg):
        raise _MustHavePROUError
    elif (deg, mod) not in _BRV_ROOTS_AND_INV_ROOTS_CACHE:
        root: int = _find_prou(mod=mod, deg=deg)
        inv_root: int = pow(base=root, exp=mod-2, mod=mod)
        powers: list[int] = [pow(base=root, exp=i, mod=mod) for i in range(deg)]
        inv_powers: list[int] = [pow(base=inv_root, exp=i, mod=mod) for i in range(deg)]
        brv_powers: list[int] = _brv_copy(val=powers)
        brv_inv_powers: list[int] = _brv_copy(val=inv_powers)
        _BRV_ROOTS_AND_INV_ROOTS_CACHE[(deg, mod)] = (brv_powers, brv_inv_powers)
    return _BRV_ROOTS_AND_INV_ROOTS_CACHE[(deg, mod)]


def _make_halfmod_logmod(mod: int) -> tuple[int, int]:
    """
    Private, cached function for make a tuple containing mod//2 and mod.bit_length().
    :param mod:
    :type mod:
    :return: A tuple containing mod//2 and mod.bit_length().
    :rtype: tuple[int, int]
    :raises _MustBeIntError: If mod is not an integer.
    :raises _MustBeOddPrimeError: If mod is not an odd prime.
    """
    if not isinstance(mod, int):
        raise _MustBeIntError
    elif not _is_odd_prime(val=mod):
        raise _MustBeOddPrimeError
    elif mod not in _MOD_HALFMOD_LOGMOD_CACHE:
        _MOD_HALFMOD_LOGMOD_CACHE[mod] = (mod//2, mod.bit_length())
    return _MOD_HALFMOD_LOGMOD_CACHE[mod]


def _cent(val: int, mod: int) -> int:
    """
    Public function to centrally reduce a value modulo mod (where mod is an odd prime) in constant-time such that the
    output satisfies -(mod//2) <= val <= mod//2.

    :param val: An integer to be reduced modulo mod.
    :type val: int
    :param mod: A modulus to reduce val by.
    :type mod: int
    :return: The value of val reduced modulo mod.
    :rtype: int
    :raises _MustBeIntError: If val or mod are not integers.
    :raises _MustBeOddPrimeError: If mod is not an odd prime.
    """
    if not isinstance(val, int) or not isinstance(mod, int):
        raise _MustBeIntError
    elif not _is_odd_prime(val=mod):
        raise _MustBeOddPrimeError
    halfmod: int
    logmod: int
    halfmod, logmod = _make_halfmod_logmod(mod=mod)
    y: int = val % mod
    w: int = y - halfmod - 1
    z: int = y - (1 + (w >> logmod)) * mod
    return z


def _cooley_tukey_ntt(val: list[int], mod: int, brv_powers: list[int]) -> list[int]:
    """
    Private in-place NTT transform of val. Input val (a list of integers in usual ordering), a modulus, and a list of
    powers of the primitive root of unity (with order 2*len(val)) in bit-reversed order.  Output the NTT transform of
    val in bit-reversed order.  See 'Speeding up the Number Theoretic Transform for Faster Ideal Lattice-Based
    Cryptography' by Longa and Naehrig (https://eprint.iacr.org/2016/504.pdf).

    :param val: A list of integers in usual ordering.
    :type val: list[int]
    :param mod: A modulus to reduce val by.
    :type mod: int
    :param brv_powers: A list of powers of the primitive root of unity in bit-reversed order.
    :type brv_powers: list[int]
    :return: The NTT transform of val in bit-reversed order.
    :rtype: list[int]
    :raises _InvalidNTTInputError: If val,mod is not ntt-friendly.
    """
    if not _is_ntt_valid(val=val, mod=mod, inv_flag=False):
        raise _InvalidNTTInputError
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
                val[j] = _cent(val=u + v, mod=mod)
                val[j + t] = _cent(val=u - v, mod=mod)
        m *= 2
    return val


def _gentleman_sande_intt(val: list[int], mod: int, brv_powers: list[int]) -> list[int]:
    """
    Private in-place INTT transform of val. Input val (a list of integers in usual ordering), a modulus, and a list of
    powers of the inverse of the primitive root of unity (with order 2*len(val)) in bit-reversed order.  Output the INTT
    transform of val in standard order.  See 'Speeding up the Number Theoretic Transform for Faster Ideal Lattice-Based
    Cryptography' by Longa and Naehrig (https://eprint.iacr.org/2016/504.pdf).

    :param val: A list of integers in usual ordering.
    :type val: list[int]
    :param mod: A modulus to reduce val by.
    :type mod: int
    :param brv_powers: A list of powers of the primitive root of unity in bit-reversed order.
    :type brv_powers: list[int]
    :return: The NTT transform of val in bit-reversed order.
    :rtype: list[int]
    :raises _InvalidNTTInputError: If val,mod is not ntt-friendly.
    """
    if not _is_ntt_valid(val=val, mod=mod, inv_flag=True):
        raise _InvalidNTTInputError
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
                val[j] = _cent(val=u + v, mod=mod)
                val[j + t] = _cent(val=(u - v) * s, mod=mod)
            j_one += 2 * t
        t *= 2
        m //= 2
    for j in range(n):
        val[j] = _cent(val=val[j] * n_inv, mod=mod)
    return val


def _ntt_poly_mult(
        f: list[int], g: list[int], mod: int, brv_powers: list[int], brv_inv_root_powers: list[int]) -> list[int]:
    """
    Private function for multiplying two coefficient representations of polynomials by:
         1. Compute the NTT transforms
         2. Compute their Hadamard product.
         3. Compute the inverse NTT transform.
    Note that since _cooley_tukey_ntt and _gentleman_sande_intt are both in-place transformations, we change f and g
    by calling this function.

    :param f: Coefficient representation of a polynomial (a list of integers).
    :type f: list[int]
    :param g: Coefficient representation of a polynomial (a list of integers).
    :type g: list[int]
    :param mod: A modulus to reduce val by.
    :type mod: int
    :param brv_powers: A list of powers of the primitive root of unity in bit-reversed order.
    :type brv_powers: list[int]
    :param brv_inv_root_powers: A list of powers of the inverse of the primitive root of unity in bit-reversed order.
    :type brv_inv_root_powers: list[int]
    :return: The coefficient representation of the product of two polynomials.
    :rtype: list[int]
    :raises _InvalidNTTInputError: If f,g,mod is not ntt-friendly.
    """
    if not _is_ntt_poly_mult_valid(f=f, g=g, mod=mod, brv_powers=brv_powers):
        raise _InvalidNTTInputError
    _cooley_tukey_ntt(val=f, mod=mod, brv_powers=brv_powers)
    _cooley_tukey_ntt(val=g, mod=mod, brv_powers=brv_powers)
    fg: list[int] = [_cent(val=x * y, mod=mod) for x, y in zip(f, g)]
    _gentleman_sande_intt(val=fg, mod=mod, brv_powers=brv_inv_root_powers)
    # Uncomment the following two lines (or make a deep copy before calling _ntt_poly_mult) to preserve f and g
    _gentleman_sande_intt(val=f, mod=mod, brv_powers=brv_inv_root_powers)
    _gentleman_sande_intt(val=g, mod=mod, brv_powers=brv_inv_root_powers)
    return fg


def _ntt(val: list[int], mod: int, inv_flag: bool) -> list[int]:
    """
    Public function for computing the Number Theoretic Transform (NTT) or its inverse in constant time.
    :param val: The list to be transformed
    :type val: list[int]
    :param mod: the modulus
    :type mod: int
    :param inv_flag: Flag indicating whether the transform is the forward transform or the inverse
    :type inv_flag: bool
    :return: the transformed list
    :rtype: list[int]
    :raises _MustBeListError: If val is not a list
    :raises _MustBeIntError: If mod is not an integer
    :raises _MustBeBoolError: If inv_flag is not a boolean
    """
    if not isinstance(val, list):
        raise _MustBeListError
    elif not isinstance(mod, int) or not all(isinstance(i, int) for i in val):
        raise _MustBeIntError
    elif not isinstance(inv_flag, bool):
        raise _MustBeBoolError
    deg: int = len(val)
    brv_powers: list[int]
    brv_inv_powers: list[int]
    brv_powers, brv_inv_powers = _brv_root_and_inv_root_powers(deg=deg, mod=mod)
    deepcopy_val = _deepcopy(val)
    if not inv_flag:
        return _cooley_tukey_ntt(val=deepcopy_val, mod=mod, brv_powers=brv_powers)
    return _gentleman_sande_intt(val=deepcopy_val, mod=mod, brv_powers=brv_inv_powers)

