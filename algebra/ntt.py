"""
algebra/ntt.py

A module for handling the Number Theoretic Transform (NTT) and its inverse in constant time.

Author: Brandon G Goodell
"""
from math import ceil as _ceil
from copy import deepcopy as _deepcopy
from typing import Dict as _Dict, Tuple as _Tuple, List as _List, Union as _Union
from algebra.errors import (_MUST_BE_INT_ERR, _MUST_BE_LIST_ERR, _MUST_BE_POS_INT_ERR,
                            _NO_PROU_FOUND_ERR, _MUST_BE_ODD_PRIME_ERR, _INVALID_NTT_INPUT_ERR,
                            _MUST_BE_POS_INT_POW_2_ERR)

# Caching dictionaries
_ODD_PRIME_CACHE: _Dict[int, bool] = {}
_HAS_PROU_CACHE: _Dict[_Tuple[int, int], bool] = {}
_POW_TWO_CACHE: _Dict[int, bool] = {}
_ROU_CACHE: _Dict[_Tuple[int, int, int], bool] = {}
_PROU_CACHE: _Dict[_Tuple[int, int, int], bool] = {}
_ROOT_INVERSE_CACHE: _Dict[_Tuple[int, int, int], bool] = {}
_FIND_PROU_CACHE: _Dict[_Tuple[int, int], _Union[int, None]] = {}
_REV_IDX_CACHE = {}
_MOD_HALFMOD_LOGMOD_CACHE: _Dict[int, _Tuple[int, int]] = {}
_BRV_ROOTS_AND_INV_ROOTS_CACHE: _Dict[_Tuple[int, int], _Tuple[_List[int], _List[int]]] = {}


# Cached boolean functions
def _is_odd_prime(val: int) -> bool:
    """
    Private, cached function for checking if x is an odd prime number.

    :param val: An integer to check for primality.
    :type val: int
    :return: True if val is an odd prime, False otherwise.
    :rtype: bool
    """
    if val not in _ODD_PRIME_CACHE:
        _ODD_PRIME_CACHE[val] = isinstance(val, int) and val >= 3
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
    if (root, inv_root, mod) not in _ROOT_INVERSE_CACHE:
        _ROOT_INVERSE_CACHE[(root, inv_root, mod)] = isinstance(root, int) and isinstance(inv_root, int)
        _ROOT_INVERSE_CACHE[(root, inv_root, mod)] = _ROOT_INVERSE_CACHE[(root, inv_root, mod)] and _is_odd_prime(val=mod)
        _ROOT_INVERSE_CACHE[(root, inv_root, mod)] = _ROOT_INVERSE_CACHE[(root, inv_root, mod)] and (root * inv_root) % mod == 1
    return _ROOT_INVERSE_CACHE[(root, inv_root, mod)]


def _is_ntt_friendly(mod: int, deg: int) -> bool:
    """
    Private, cached function foir checking if a modulus has a primitive root of unity of order 2*deg.

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


# Functionality
# find_prou - find a primitive root of unity
_GENERAL_FIND_PROU_ERR: str = "find prou validation error"
_NOT_NTT_FRIENDLY_ERR: str = "degree-modulus pair not ntt-friendly"


class FindProuError(Exception):
    def __init__(self, message: str = _GENERAL_FIND_PROU_ERR) -> None:
        self.message = message
        super().__init__(message)


class FindProuMustBeNTTFriendlyError(FindProuError):
    def __init__(self, message: str = _NOT_NTT_FRIENDLY_ERR) -> None:
        self.message = message
        super().__init__(message)


class FindProuNoProuFoundError(FindProuError):
    def __init__(self, message: str = _NO_PROU_FOUND_ERR) -> None:
        self.message = message
        super().__init__(message)


def _check_find_prou(mod: int, deg: int):
    """
    Raise MustHaveProuError if mod and deg are not ntt-friendly.

    :param mod: A modulus to check for having a primitive root of unity of order 2*deg.
    :type mod: int
    :param deg: The degree of the polynomial.
    :type deg: int
    :raises MustHaveProuError: If mod and deg are not ntt-friendly.
    """
    if not _is_ntt_friendly(mod=mod, deg=deg):
        raise FindProuMustBeNTTFriendlyError


def _find_prou(mod: int, deg: int) -> int:
    """
    Private, cached function to compute primitive root of unity of order 2*deg modulo mod where mod and deg are
    ntt-friendly. Computes with a naive loop that first checks 2, then 3, then 4, and so on.

    :param mod: A modulus to check for having a primitive root of unity of order 2*deg.
    :type mod: int
    :param deg: The degree of the polynomial.
    :type deg: int
    :return: A primitive root of unity of order 2*deg modulo mod.
    :rtype: int
    :raises RuntimeError: If no primitive root is found despite that mod and deg are ntt-friendly (should not happen).
    """
    if (mod, deg) not in _FIND_PROU_CACHE:
        _FIND_PROU_CACHE[(mod, deg)] = None
        r: int = 2
        while r < mod and not _is_prou(root=r, mod=mod, deg=deg):
            r += 1
        if r is None or not _is_prou(root=r, mod=mod, deg=deg):
            raise FindProuNoProuFoundError
        _FIND_PROU_CACHE[(mod, deg)] = r
    return _FIND_PROU_CACHE[(mod, deg)]


def find_prou(mod: int, deg: int) -> int:
    """
    Public function to compute primitive root of unity of order 2*deg modulo mod where mod and deg are ntt-friendly.

    :param mod: A modulus to check for having a primitive root of unity of order 2*deg.
    :type mod: int
    :param deg: The degree of the polynomial.
    :type deg: int
    :return: A primitive root of unity of order 2*deg modulo mod.
    :rtype: int
    :raises MustHaveProuError: If mod and deg are not ntt-friendly.
    """
    _check_find_prou(mod=mod, deg=deg)
    return _find_prou(mod=mod, deg=deg)


# _brv_indices - generate the permutation to bit-reverse indices in a list of power-2 length
_GENERAL_BRV_INDICES_ERR : str = "BRV indices validation error"


class BRVIndicesValidationError(Exception):
    def __init__(self, message: str = _GENERAL_BRV_INDICES_ERR) -> None:
        self.message = message
        super().__init__(message)


class BRVIndicesMustBeIntError(BRVIndicesValidationError):
    def __init__(self, message: str = _MUST_BE_INT_ERR) -> None:
        self.message = message
        super().__init__(message)


class BRVIndicesMustBeStrictlyPostiive(BRVIndicesValidationError):
    def __init__(self, message: str = _MUST_BE_POS_INT_ERR) -> None:
        self.message = message
        super().__init__(message)


class BRVIndicesMustBePowerOfTwoError(BRVIndicesValidationError):
    def __init__(self, message: str = _MUST_BE_POS_INT_POW_2_ERR) -> None:
        self.message = message
        super().__init__(message)


def _is_brv_indices_valid(val: int) -> bool:
    """
    Private function for checking if input to bit-reversed indices is valid.

    :param val:
    :type val: int
    :return: boolean indicating whether input val is a positive power of 2.
    :rtype: bool
    """
    return isinstance(val, int) and _is_pos_pow_two(val=val)


def _check_brv_indices(val: int):
    """
    Private function for checking if input to bit-reversed indices is valid.

    :param val:
    :type val: int
    :raises BRVIndicesValidationError: If val is not a positive power of 2.
    """
    if not isinstance(val, int):
        raise BRVIndicesMustBeIntError
    elif val < 1:
        raise BRVIndicesMustBeStrictlyPostiive
    elif not _is_pos_pow_two(val=val):
        raise BRVIndicesMustBePowerOfTwoError


def __brv_indices(val: int) -> list[int]:
    """
    Private, cached function for computing bit-reversed indices of a list of length val.

    :param val:
    :type val: int
    :return: list of bit-reversal-permuted indices.
    """
    if val not in _REV_IDX_CACHE:
        k: int = val.bit_length() - 1
        _REV_IDX_CACHE[val] = [int(bin(i)[2:].zfill(k)[::-1], 2) for i in range(val)]
    return _REV_IDX_CACHE[val]


def _brv_indices(val: int) -> list[int]:
    """
    Private, cached function for computing bit-reversed indices of a list of length val.

    :param val: The length of the list.
    :type val: int
    :return: A list of bit-reversed indices.
    :rtype: List[int]
    :raises BRVIndicesValidationError: If val is not a positive power of 2.
    """
    _check_brv_indices(val=val)
    return __brv_indices(val=val)


# brv_copy - compute a copy of a list, permuting indices by bit-reversal
_GENERAL_BRV_COPY_ERR: str = "BRV copy validation error"


class BRVCopyValidationError(Exception):
    def __init__(self, message: str = _GENERAL_BRV_COPY_ERR) -> None:
        self.message = message
        super().__init__(message)


class BRVCopyMustBeListError(BRVCopyValidationError):
    def __init__(self, message: str = _MUST_BE_LIST_ERR) -> None:
        self.message = message
        super().__init__(message)


def _check_brv_copy(val: list):
    if not isinstance(val, list):
        raise BRVCopyMustBeListError


def _brv_copy(val: list) -> list:
    """
    Private function for computing a copy of 'val', permuting indices by bit-reversal. Do not cache.

    :param val: A list of integers.
    :type val: list
    :return: A list of integers permuted by bit-reversal.
    :rtype: list
    :raises TypeError: If val is not a list.
    """
    brvd_indices: _List[int] = _brv_indices(val=len(val))
    return [val[i] for i in brvd_indices]


def brv_copy(val: list) -> list:
    """
    Public function for computing a copy of 'val', permuting indices by bit-reversal. Do not cache.

    :param val: A list of integers.
    :type val: list
    :return: A list of integers permuted by bit-reversal.
    :rtype: list
    :raises TypeError: If val is not a list.
    """
    _check_brv_copy(val=val)
    return _brv_copy(val=val)


# _brv_root_and_inv_root_powers
_GENERAL_BRV_ROOT_AND_INV_ROOT_POWERS_ERR: str = "BRV root and inv root powers validation error"


class BRVRootAndInvRootPowers(Exception):
    def __init__(self, message: str = _GENERAL_BRV_ROOT_AND_INV_ROOT_POWERS_ERR) -> None:
        self.message = message
        super().__init__(message)


class BRVRootAndInvRootMustBeNTTFriendlyError(BRVRootAndInvRootPowers):
    def __init__(self, message: str = _NOT_NTT_FRIENDLY_ERR) -> None:
        self.message = message
        super().__init__(message)


def _check_brv_root_and_inv_root_powers(deg: int, mod: int):
    if not _is_ntt_friendly(mod=mod, deg=deg):
        raise BRVRootAndInvRootMustBeNTTFriendlyError


def __brv_root_and_inv_root_powers(deg: int, mod: int) -> _Tuple[_List[int], _List[int]]:
    if (deg, mod) not in _BRV_ROOTS_AND_INV_ROOTS_CACHE:
        root: int = find_prou(mod=mod, deg=deg)
        inv_root: int = pow(base=root, exp=mod - 2, mod=mod)
        powers: _List[int] = [pow(base=root, exp=i, mod=mod) for i in range(deg)]
        inv_powers: _List[int] = [pow(base=inv_root, exp=i, mod=mod) for i in range(deg)]
        brv_powers: _List[int] = brv_copy(val=powers)
        brv_inv_powers: _List[int] = brv_copy(val=inv_powers)
        _BRV_ROOTS_AND_INV_ROOTS_CACHE[(deg, mod)] = (brv_powers, brv_inv_powers)
    return _BRV_ROOTS_AND_INV_ROOTS_CACHE[(deg, mod)]


def _brv_root_and_inv_root_powers(deg: int, mod: int) -> _Tuple[_List[int], _List[int]]:
    """
    Private, cached function for computing the primitive root of unity of order 2*deg modulo mod and its inverse, and
    returning their powers.

    :param deg: The degree of the polynomial.
    :type deg: int
    :param mod: A modulus to check for having a primitive root of unity of order 2*deg.
    :type mod: int
    :return: A tuple containing the powers of the primitive root of unity and its inverse.
    :rtype: Tuple[List[int], List[int]]
    :raises ValueError: If mod and deg are not ntt-friendly.
    """
    _check_brv_root_and_inv_root_powers(deg=deg, mod=mod)
    return __brv_root_and_inv_root_powers(deg=deg, mod=mod)


# _make_halfmod_logmod
_GENERAL_MAKE_HALFMOD_LOGMOD_ERR: str = "make_halfmod_logmod validation error"


class MakeHalfmodLogmodError(Exception):
    def __init__(self, message: str = _GENERAL_MAKE_HALFMOD_LOGMOD_ERR) -> None:
        self.message = message
        super().__init__(message)


class MakeHalfmodLogmodMustBePosIntError(MakeHalfmodLogmodError):
    def __init__(self, message: str = _MUST_BE_POS_INT_ERR) -> None:
        self.message = message
        super().__init__(message)


def _check_make_halfmod_logmod(mod: int):
    if not isinstance(mod, int) or (isinstance(mod, int) and mod < 1):
        raise MakeHalfmodLogmodMustBePosIntError


def _make_halfmod_logmod(mod: int) -> _Tuple[int, int]:
    """
    Private, cached function for make a tuple containing mod//2 and mod.bit_length().
    :param mod:
    :type mod:
    :return: A tuple containing mod//2 and mod.bit_length().
    :rtype: Tuple[int, int]
    """
    if mod not in _MOD_HALFMOD_LOGMOD_CACHE:
        _MOD_HALFMOD_LOGMOD_CACHE[mod] = (mod//2, mod.bit_length())
    return _MOD_HALFMOD_LOGMOD_CACHE[mod]


def make_halfmod_logmod(mod: int) -> _Tuple[int, int]:
    """
    Public function for make a tuple containing mod//2 and mod.bit_length().
    :param mod:
    :type mod:
    :return: A tuple containing mod//2 and mod.bit_length().
    :rtype: Tuple[int, int]
    :raises ValueError: If mod is not a positive integer.
    """
    _check_make_halfmod_logmod(mod=mod)
    return _make_halfmod_logmod(mod=mod)


# cent
_GENERAL_CENT_ERR: str = "cent validation error"


class CentError(Exception):
    def __init__(self, message: str = _GENERAL_CENT_ERR) -> None:
        self.message = message
        super().__init__(message)


class CentMustBeIntError(CentError):
    def __init__(self, message: str = _MUST_BE_INT_ERR) -> None:
        self.message = message
        super().__init__(message)


class CentMustBeOddPrimeError(CentError):
    def __init__(self, message: str = _MUST_BE_ODD_PRIME_ERR) -> None:
        self.message = message
        super().__init__(message)


def _check_cent(val: int, mod: int):
    if not isinstance(val, int):
        raise CentMustBeIntError
    elif not _is_odd_prime(val=mod):
        raise CentMustBeOddPrimeError


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
    :raises TypeError: If val is not an integer.
    :raises ValueError: If mod is not an odd prime.
    """
    halfmod: int
    logmod: int
    halfmod, logmod = _make_halfmod_logmod(mod=mod)
    y: int = val % mod
    w: int = y - halfmod - 1
    z: int = y - (1 + (w >> logmod)) * mod
    return z


def cent(val: int, mod: int) -> int:
    """
    Public function to centrally reduce a value modulo mod (where mod is an odd prime) in constant-time such that the
    output satisfies -(mod//2) <= val <= mod//2.

    :param val: An integer to be reduced modulo mod.
    :type val: int
    :param mod: A modulus to reduce val by.
    :type mod: int
    :return: The value of val reduced modulo mod.
    :rtype: int
    :raises TypeError: If val is not an integer.
    :raises ValueError: If mod is not an odd prime.
    """
    _check_cent(val=val, mod=mod)
    return _cent(val=val, mod=mod)


#
_GENERAL_NTT_ERROR: str = "ntt validation error"


class NTTError(Exception):
    def __init__(self, message: str = _GENERAL_NTT_ERROR) -> None:
        self.message = message
        super().__init__(message)


class NTTMustBeListError(NTTError):
    def __init__(self, message: str = _MUST_BE_LIST_ERR) -> None:
        self.message = message
        super().__init__(message)


class NTTMustBeIntError(NTTError):
    def __init__(self, message: str = _MUST_BE_INT_ERR) -> None:
        self.message = message
        super().__init__(message)


class NTTMustBeBoolError(NTTError):
    def __init__(self, message: str = _MUST_BE_INT_ERR) -> None:
        self.message = message
        super().__init__(message)


class NTTMustBeNTTFriendlyError(NTTError):
    def __init__(self, message: str = _NOT_NTT_FRIENDLY_ERR) -> None:
        self.message = message
        super().__init__(message)


def _check_ntt(val: list[int], mod: int, inv_flag: bool):
    if not isinstance(val, list):
        raise NTTMustBeListError
    elif not all(isinstance(x, int) for x in val+[mod]):
        raise NTTMustBeIntError
    elif not isinstance(inv_flag, bool):
        raise NTTMustBeBoolError
    elif not _is_ntt_friendly(mod=mod, deg=len(val)):
        raise NTTMustBeNTTFriendlyError


def _cooley_tukey_ntt(val: _List[int], mod: int, brv_powers: _List[int]) -> _List[int]:
    """
    Private in-place NTT transform of val. Input val (a list of integers in usual ordering), a modulus, and a list of powers of
    the primitive root of unity (with order 2*len(val)) in bit-reversed order.  Output the NTT transform of val in
    bit-reversed order.  See 'Speeding up the Number Theoretic Transform for Faster Ideal Lattice-Based Cryptography' by
    Longa and Naehrig (https://eprint.iacr.org/2016/504.pdf).

    :param val: A list of integers in usual ordering.
    :type val: List[int]
    :param mod: A modulus to reduce val by.
    :type mod: int
    :param brv_powers: A list of powers of the primitive root of unity in bit-reversed order.
    :type brv_powers: List[int]
    :return: The NTT transform of val in bit-reversed order.
    :rtype: List[int]
    :raises TypeError: If val is not ntt-friendly.
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
            s: int = brv_powers[m + i]
            for j in range(j_one, j_two + 1):
                u, v = val[j], val[j + t] * s
                val[j] = cent(val=u + v, mod=mod)
                val[j + t] = cent(val=u - v, mod=mod)
        m *= 2
    return val


def _gentleman_sande_intt(val: _List[int], mod: int, brv_powers: _List[int]) -> _List[int]:
    """
    Private in-place INTT transform of val. Input val (a list of integers in usual ordering), a modulus, and a list of powers of
    the inverse of the primitive root of unity (with order 2*len(val)) in bit-reversed order.  Output the INTT transform
    of val in standard order.  See 'Speeding up the Number Theoretic Transform for Faster Ideal Lattice-Based
    Cryptography' by Longa and Naehrig (https://eprint.iacr.org/2016/504.pdf).

    :param val: A list of integers in usual ordering.
    :type val: List[int]
    :param mod: A modulus to reduce val by.
    :type mod: int
    :param brv_powers: A list of powers of the primitive root of unity in bit-reversed order.
    :type brv_powers: List[int]
    :return: The NTT transform of val in bit-reversed order.
    :rtype: List[int]
    :raises TypeError: If val is not ntt-friendly.
    """
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
                val[j] = cent(val=u + v, mod=mod)
                val[j + t] = cent(val=(u - v) * s, mod=mod)
            j_one += 2 * t
        t *= 2
        m //= 2
    for j in range(n):
        val[j] = cent(val=val[j] * n_inv, mod=mod)
    return val


def _ntt(val: _List[int], mod: int, inv_flag: bool) -> _List[int]:
    """
    Public wrapper for NTT/INTT transforms.

    :param val: A list of integers.
    :type val: List[int]
    :param mod: The modulus.
    :type mod: int
    :param inv_flag: A flag indicating forward or inverse.
    :type inv_flag: bool
    :return: The NTT or INTT transform of val.
    :rtype: List[int]
    """
    deg: int = len(val)
    brv_powers: _List[int]
    brv_inv_powers: _List[int]
    brv_powers, brv_inv_powers = _brv_root_and_inv_root_powers(deg=deg, mod=mod)
    deepcopy_val = _deepcopy(val)
    if not inv_flag:
        return _cooley_tukey_ntt(val=deepcopy_val, mod=mod, brv_powers=brv_powers)
    return _gentleman_sande_intt(val=deepcopy_val, mod=mod, brv_powers=brv_inv_powers)


def ntt(val: list[int], mod: int, inv_flag: bool) -> list[int]:
    _check_ntt(val=val, mod=mod, inv_flag=inv_flag)
    return _ntt(val=val, mod=mod, inv_flag=inv_flag)





#
# def _is_ntt_poly_mult_valid(f: _List[int], g: _List[int], mod: int, brv_powers: _List[int]):
#     """
#     Private function for checking if inputs to NTT polynomial multiplication are valid.
#
#     :param f: Coefficient representation of a polynomial (a list of integers).
#     :type f: List[int]
#     :param g: Coefficient representation of a polynomial (a list of integers).
#     :type g: List[int]
#     :param mod: A modulus to reduce val by.
#     :type mod: int
#     :param brv_powers: A list of powers of the primitive root of unity in bit-reversed order.
#     :type brv_powers: List[int]
#     :return: True if the inputs to NTT polynomial multiplication are valid, False otherwise.
#     :rtype: bool
#     """
#     deg: int = len(f)
#     f_is_list = isinstance(f, list)
#     g_is_list = isinstance(g, list)
#     mod_is_int = isinstance(mod, int)
#     brv_powers_is_list = isinstance(brv_powers, list)
#     f_and_g_are_ints = all(isinstance(i, int) for i in f+g)
#     f_and_g_have_len_deg = len(f) == len(g) == deg
#     have_prou = _is_ntt_friendly(mod=mod, deg=deg)
#     powers: _List[int] = brv_copy(val=brv_powers)
#     root: int = powers[1]
#     inv_root: int = pow(base=root, exp=mod - 2, mod=mod)
#     have_correct_root = _is_prou(root=root, mod=mod, deg=deg)
#     have_correct_inverse_root = _is_root_inverse(root=root, inv_root=inv_root, mod=mod)
#     ntt_is_valid_for_f = _is_ntt_valid(val=f, mod=mod, inv_flag=False)
#     ntt_is_valid_for_g = _is_ntt_valid(val=g, mod=mod, inv_flag=False)
#     return all([f_is_list, g_is_list, mod_is_int, brv_powers_is_list, f_and_g_are_ints, f_and_g_have_len_deg, have_prou,
#                 have_correct_root, have_correct_inverse_root, ntt_is_valid_for_f, ntt_is_valid_for_g])
#
#
# def _ntt_poly_mult(f: _List[int], g: _List[int], mod: int, brv_powers: _List[int], brv_inv_root_powers: _List[int]) -> _List[int]:
#     """
#     Private function for multiplying two coefficient representations of polynomials by:
#          1. Compute the NTT transforms
#          2. Compute their Hadamard product.
#          3. Compute the inverse NTT transform.
#     Note that since _cooley_tukey_ntt and _gentleman_sande_intt are both in-place transformations, we change f and g
#     by calling this function.
#
#     :param f: Coefficient representation of a polynomial (a list of integers).
#     :type f: List[int]
#     :param g: Coefficient representation of a polynomial (a list of integers).
#     :type g: List[int]
#     :param mod: A modulus to reduce val by.
#     :type mod: int
#     :param brv_powers: A list of powers of the primitive root of unity in bit-reversed order.
#     :type brv_powers: List[int]
#     :param brv_inv_root_powers: A list of powers of the inverse of the primitive root of unity in bit-reversed order.
#     :type brv_inv_root_powers: List[int]
#     :return: The coefficient representation of the product of two polynomials.
#     :rtype: List[int]
#     :raises TypeError: If _is_ntt_poly_mult_valid returns False.
#     """
#     if not _is_ntt_poly_mult_valid(f=f, g=g, mod=mod, brv_powers=brv_powers):
#         raise TypeError(_INVALID_NTT_INPUT_ERR)
#     _cooley_tukey_ntt(val=f, mod=mod, brv_powers=brv_powers)
#     _cooley_tukey_ntt(val=g, mod=mod, brv_powers=brv_powers)
#     fg: _List[int] = [cent(val=x * y, mod=mod) for x, y in zip(f, g)]
#     _gentleman_sande_intt(val=fg, mod=mod, brv_powers=brv_inv_root_powers)
#     # Uncomment the following two lines (or make a deep copy before calling _ntt_poly_mult) to preserve f and g
#     _gentleman_sande_intt(val=f, mod=mod, brv_powers=brv_inv_root_powers)
#     _gentleman_sande_intt(val=g, mod=mod, brv_powers=brv_inv_root_powers)
#     return fg
#
#
# def ntt_poly_mult(f: _List[int], g: _List[int], mod: int) -> _List[int]:
#     """
#     Public wrapper of _ntt_poly_mult which infers deg from len(f), computes the powers of the primitive root of unity
#     and its inverse, and then calls _ntt_poly_mult.
#
#     :param f: Coefficient representation of a polynomial (a list of integers).
#     :type f: List[int]
#     :param g: Coefficient representation of a polynomial (a list of integers).
#     :type g: List[int]
#     :param mod: A modulus to reduce val by.
#     :type mod: int
#     :return: The coefficient representation of the product of two polynomials.
#     :rtype: List[int]
#     """
#     deg: int = len(f)
#     brv_powers: _List[int]
#     brv_inv_powers: _List[int]
#     brv_powers, brv_inv_powers = _brv_root_and_inv_root_powers(deg=deg, mod=mod)
#     return _ntt_poly_mult(f=f, g=g, mod=mod, brv_powers=brv_powers, brv_inv_root_powers=brv_inv_powers)
#
