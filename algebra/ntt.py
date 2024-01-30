"""
The ntt module handles the Number Theoretic Transform (NTT) and its inverse in constant time.
"""
from math import ceil as _ceil
from typing import Dict, Tuple, List, Union
from algebra.errors import (_MUST_BE_INT_ERR, _MUST_BE_LIST_ERR, _MUST_BE_POS_INT_ERR, _MUST_HAVE_PROU_ERR, 
                            _NO_PROU_FOUND_ERR, _MUST_BE_ODD_PRIME_ERR, _INVALID_NTT_INPUT_ERR)

# Caching dictionaries
_ODD_PRIME_CACHE: Dict[int, bool] = {}
_HAS_PROU_CACHE: Dict[Tuple[int, int], bool] = {}
_POW_TWO_CACHE: Dict[int, bool] = {}
_ROU_CACHE: Dict[Tuple[int, int, int], bool] = {}
_PROU_CACHE: Dict[Tuple[int, int, int], bool] = {}
_FIND_PROU_CACHE: Dict[Tuple[int, int], Union[int, None]] = {}
_REV_IDX_CACHE = {}
_MOD_HALFMOD_LOGMOD_CACHE: Dict[int, Tuple[int, int, int]] = {}
_BRV_ROOTS_AND_INV_ROOTS_CACHE: Dict[Tuple[int, int], Tuple[List[int], List[int]]] = {}


# Functions
def _is_odd_prime(val: int) -> bool:
    """
    Check if x is an odd prime number.
    """
    if val not in _ODD_PRIME_CACHE:
        _ODD_PRIME_CACHE[val] = isinstance(val, int)
        _ODD_PRIME_CACHE[val] = _ODD_PRIME_CACHE[val] and val >= 3
        _ODD_PRIME_CACHE[val] = _ODD_PRIME_CACHE[val] and all(val % i != 0 for i in range(2, _ceil(val ** 0.5) + 1))
    return _ODD_PRIME_CACHE[val]


def _is_pos_pow_two(val: int) -> bool:
    """
    Check if a number is a power of 2.
    """
    if val not in _POW_TWO_CACHE:
        _POW_TWO_CACHE[val] = isinstance(val, int)
        _POW_TWO_CACHE[val] = _POW_TWO_CACHE[val] and val >= 1
        _POW_TWO_CACHE[val] = _POW_TWO_CACHE[val] and val & (val - 1) == 0
    return _POW_TWO_CACHE[val]


def _is_ntt_friendly(mod: int, deg: int) -> bool:
    """
    Check if a modulus has a primitive root of unity of order 2*deg
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
    Check if val is a root of unity of order 2*deg modulo mod for ntt-friendly mod and deg
    """
    if (root, mod, deg) not in _ROU_CACHE:
        _ROU_CACHE[(root, mod, deg)] = _is_ntt_friendly(mod=mod, deg=deg)
        _ROU_CACHE[(root, mod, deg)] = _ROU_CACHE[(root, mod, deg)] and isinstance(root, int)
        _ROU_CACHE[(root, mod, deg)] = _ROU_CACHE[(root, mod, deg)] and pow(base=root, exp=2 * deg, mod=mod) == 1
    return _ROU_CACHE[(root, mod, deg)]


def _is_prou(root: int, mod: int, deg: int) -> bool:
    """
    Check if val is a primitive root of order 2*deg modulo modulus.
    """
    if (root, mod, deg) not in _PROU_CACHE:
        _PROU_CACHE[(root, mod, deg)] = _is_rou(root=root, mod=mod, deg=deg)
        _PROU_CACHE[(root, mod, deg)] = _PROU_CACHE[(root, mod, deg)] and all(pow(base=root, exp=i, mod=mod) != 1 for i in range(1, 2*deg))
    return _PROU_CACHE[(root, mod, deg)]


def _is_root_inverse(root: int, inv_root: int, mod: int) -> bool:
    """
    Check if root and inv_root are multiplicative inverses modulo mod for an odd prime mod
    """
    return isinstance(root, int) and isinstance(inv_root, int) and _is_odd_prime(val=mod) and (root * inv_root) % mod == 1


def _find_prou(mod: int, deg: int) -> int:
    """
    Find a primitive root of order 2*deg modulo mod where mod and deg are ntt-friendly.
    Naive loop that first checks 2, then 3, then 4...
    """
    if not _is_ntt_friendly(mod=mod, deg=deg):
        raise ValueError(_MUST_HAVE_PROU_ERR)
    elif (mod, deg) not in _FIND_PROU_CACHE:
        _FIND_PROU_CACHE[(mod, deg)] = None
        r: int = 2
        while r < mod and not _is_prou(root=r, mod=mod, deg=deg):
            r += 1
        if r is None or not _is_prou(root=r, mod=mod, deg=deg):
            raise RuntimeError(_NO_PROU_FOUND_ERR)
        _FIND_PROU_CACHE[(mod, deg)] = r
    return _FIND_PROU_CACHE[(mod, deg)]


def _brv_indices(val: int) -> List[int]:
    """
    Compute bit-reversed indices of a list of length val.
    """
    if not isinstance(val, int):
        raise TypeError(_MUST_BE_INT_ERR)
    elif not _is_pos_pow_two(val=val):
        raise ValueError(_MUST_BE_POS_INT_ERR)
    elif val not in _REV_IDX_CACHE:
        k: int = val.bit_length() - 1
        _REV_IDX_CACHE[val] = [int(bin(i)[2:].zfill(k)[::-1], 2) for i in range(val)]
    return _REV_IDX_CACHE[val]


def _brv_copy(val: list) -> list:
    """
    Permute indices by bit-reversal.
    """
    if not isinstance(val, list):
        raise TypeError(_MUST_BE_LIST_ERR)
    brvd_indices: List[int] = _brv_indices(val=len(val))
    return [val[i] for i in brvd_indices]


def _brv_root_and_inv_root_powers(deg: int, mod: int) -> Tuple[List[int], List[int]]:
    """
    Compute primitive root of unity of order 2*deg modulo mod, and its inverse, and then return their powers
    """
    if not _is_ntt_friendly(mod=mod, deg=deg):
        raise ValueError(_MUST_HAVE_PROU_ERR)
    elif (deg, mod) not in _BRV_ROOTS_AND_INV_ROOTS_CACHE:
        root: int = _find_prou(mod=mod, deg=deg)
        inv_root: int = pow(base=root, exp=mod-2, mod=mod)
        powers: List[int] = [pow(base=root, exp=i, mod=mod) for i in range(deg)]
        inv_powers: List[int] = [pow(base=inv_root, exp=i, mod=mod) for i in range(deg)]
        brv_powers: List[int] = _brv_copy(val=powers)
        brv_inv_powers: List[int] = _brv_copy(val=inv_powers)
        _BRV_ROOTS_AND_INV_ROOTS_CACHE[(deg, mod)] = (brv_powers, brv_inv_powers)
    return _BRV_ROOTS_AND_INV_ROOTS_CACHE[(deg, mod)]


def _cent(val: int, mod: int) -> int:
    """
    Centrally reduce a value modulo mod in constant time where mod is an odd prime
    Output val satisfies -(mod//2) <= val <= mod//2.
    """
    if not isinstance(val, int):
        raise TypeError(_MUST_BE_INT_ERR)
    elif not _is_odd_prime(val=mod):
        raise ValueError(_MUST_BE_ODD_PRIME_ERR)
    y: int = val % mod
    w: int = y - (mod//2) - 1
    z: int = y - (1 + (w >> mod.bit_length())) * mod
    return z


def _is_ntt_valid(val: List[int], mod: int, inv_flag: bool) -> bool:
    """
    Check if input to NTT is valid.
    """
    val_is_list: bool = isinstance(val, list)
    deg: int = len(val)
    deg_mod_is_ntt_friendly: bool = _is_ntt_friendly(mod=mod, deg=deg)
    inv_flag_is_bool: bool = isinstance(inv_flag, bool)
    entries_in_val_are_ints: bool = all(isinstance(i, int) for i in val)
    brv_powers: List[int]
    brv_powers, brv_inv_powers = _brv_root_and_inv_root_powers(deg=deg, mod=mod)
    powers: List[int] = _brv_copy(val=brv_powers)
    inv_powers: List[int] = _brv_copy(val=brv_inv_powers)
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


def _cooley_tukey_ntt(val: List[int], mod: int, brv_powers: List[int]) -> List[int]:
    """
    Input val, a list of n := len(val) integers in usual ordering, a modulus that is a prime such that
    (modulus-1) % (2*deg) == 0, a root_order == 2*deg, and a list of integers, root_powers, which are powers of a primitive
    root of unity with order 2*deg stored in bit-reversed order.  Output the NTT of val in bit-reversed order.

    In-place computation, iterative implementation of the Cooley-Tukey butterfly as defined in 'Speeding up the Number
    Theoretic Transform for Faster Ideal Lattice-Based Cryptography' by Longa and Naehrig
    (https://eprint.iacr.org/2016/504.pdf).
    """
    if not _is_ntt_valid(val=val, mod=mod, inv_flag=False):
        raise TypeError(_INVALID_NTT_INPUT_ERR)
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


def _gentleman_sande_intt(val: List[int], mod: int, brv_powers: List[int]) -> List[int]:
    """
    Input val, a list of n := len(val) integers in usual ordering, a modulus that is a prime such that
    (modulus-1) % (2*n) == 0, a root_order == 2*deg, and a list of integers, inv_root_powers, which are powers of the
    inverse of the primitive root of unity with order 2*deg used to compute the forward transform in bit-reversed
    order.  Output the INTT of val in standard order.

    In-place computation, iterative implementation of the Gentleman-Sande butterfly as in 'Speeding up the Number
    Theoretic Transform for Faster Ideal Lattice-Based Cryptography' by Longa and Naehrig
    (https://eprint.iacr.org/2016/504.pdf).
    """
    if not _is_ntt_valid(val=val, mod=mod, inv_flag=True):
        raise TypeError(_INVALID_NTT_INPUT_ERR)
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


def _is_ntt_poly_mult_valid(f: List[int], g: List[int], mod: int, brv_powers: List[int]):
    """
    Check if inputs to NTT polynomial multiplication are valid.
    """
    deg: int = len(f)
    f_is_list = isinstance(f, list)
    g_is_list = isinstance(g, list)
    mod_is_int = isinstance(mod, int)
    brv_powers_is_list = isinstance(brv_powers, list)
    f_and_g_are_ints = all(isinstance(i, int) for i in f+g)
    f_and_g_have_len_deg = len(f) == len(g) == deg
    have_prou = _is_ntt_friendly(mod=mod, deg=deg)
    powers: List[int] = _brv_copy(val=brv_powers)
    root: int = powers[1]
    inv_root: int = pow(base=root, exp=mod - 2, mod=mod)
    have_correct_root = _is_prou(root=root, mod=mod, deg=deg)
    have_correct_inverse_root = _is_root_inverse(root=root, inv_root=inv_root, mod=mod)
    ntt_is_valid_for_f = _is_ntt_valid(val=f, mod=mod, inv_flag=False)
    ntt_is_valid_for_g = _is_ntt_valid(val=g, mod=mod, inv_flag=False)
    return all([f_is_list, g_is_list, mod_is_int, brv_powers_is_list, f_and_g_are_ints, f_and_g_have_len_deg, have_prou, have_correct_root, have_correct_inverse_root, ntt_is_valid_for_f, ntt_is_valid_for_g])


def _ntt_poly_mult(f: List[int], g: List[int], mod: int, brv_powers: List[int], brv_inv_root_powers: List[int]) -> List[int]:
    """
    Multiply two coefficient representations of polynomials by first computing their NTTs
    and then their Hadamard product before inverting the NTT.
    """
    if not _is_ntt_poly_mult_valid(f=f, g=g, mod=mod, brv_powers=brv_powers):
        raise TypeError(_INVALID_NTT_INPUT_ERR)
    _cooley_tukey_ntt(val=f, mod=mod, brv_powers=brv_powers)
    _cooley_tukey_ntt(val=g, mod=mod, brv_powers=brv_powers)
    fg: List[int] = [_cent(val=x * y, mod=mod) for x, y in zip(f, g)]
    _gentleman_sande_intt(val=fg, mod=mod, brv_powers=brv_inv_root_powers)
    _gentleman_sande_intt(val=f, mod=mod, brv_powers=brv_inv_root_powers)
    _gentleman_sande_intt(val=g, mod=mod, brv_powers=brv_inv_root_powers)
    return fg
