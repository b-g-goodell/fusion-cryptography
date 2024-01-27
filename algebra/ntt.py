"""
The ntt module handles the Number Theoretic Transform (NTT) and its inverse in constant time.
"""
from math import ceil as _ceil
from typing import Dict as _Dict, Tuple as _Tuple, List as _List, Union as _Union
from algebra.errors import (_MUST_BE_INT_ERR, _MUST_BE_LIST_ERR, _MUST_BE_LIST_W_POW_2_LEN_ERR, _MUST_BE_POS_INT_ERR,
                            _MUST_BE_INT_GEQ_3_ERR, _MUST_BE_HALF_FLOORED_ERR, _MUST_BE_BIT_LEN_ERR, _MUST_HAVE_PROU_ERR,
                            _NO_PROU_FOUND_ERR, _MUST_BE_BOOL_ERR, _MUST_BE_ODD_PRIME_ERR, _MUST_BE_CORRECT_ROOT_ERR,
                            _MUST_BE_CORRECT_INVERSE_ROOT_ERR, _DEGREE_MISMATCH_ERR,
                            _MUST_CONSTRUCT_BRV_POWERS_CORRECTLY_ERR, _INVALID_NTT_INPUT_ERR)

# Caching dictionaries
_ODD_PRIME_CACHE: _Dict[int, bool] = {}
_HAS_PROU_CACHE: _Dict[_Tuple[int, int], bool] = {}
_POW_TWO_CACHE: _Dict[int, bool] = {}
_ROU_CACHE: _Dict[_Tuple[int, int, int], bool] = {}
_PROU_CACHE: _Dict[_Tuple[int, int, int], bool] = {}
_FIND_PROU_CACHE: _Dict[_Tuple[int, int], _Union[int, None]] = {}
_REV_IDX_CACHE = {}
_MOD_HALFMOD_LOGMOD_CACHE: _Dict[int, _Tuple[int, int, int]] = {}
_DERIVED_PARAMS_CACHE: _Dict[_Tuple[int, int, bool], _List[int]] = {}


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


def _has_prou(mod: int, deg: int) -> bool:
    """
    Check if a modulus has a primitive root of unity of a given order.
    """
    if (mod, deg) not in _HAS_PROU_CACHE:
        _HAS_PROU_CACHE[(mod, deg)] = isinstance(mod, int)
        _HAS_PROU_CACHE[(mod, deg)] = _HAS_PROU_CACHE[(mod, deg)] and isinstance(deg, int)
        _HAS_PROU_CACHE[(mod, deg)] = _HAS_PROU_CACHE[(mod, deg)] and mod >= 3
        _HAS_PROU_CACHE[(mod, deg)] = _HAS_PROU_CACHE[(mod, deg)] and deg >= 1
        _HAS_PROU_CACHE[(mod, deg)] = _HAS_PROU_CACHE[(mod, deg)] and (mod - 1) % (2 * deg) == 0
    return _HAS_PROU_CACHE[(mod, deg)]


def _is_pow_two_geq_two(val: int) -> bool:
    """
    Check if a number is a power of two.
    """
    if val not in _POW_TWO_CACHE:
        _POW_TWO_CACHE[val] = isinstance(val, int) and val >= 2 and val & (val - 1) == 0
    return _POW_TWO_CACHE[val]


def _bit_reverse_len_to_idxs(val: int) -> _List[int]:
    """
    Compute bit-reversed indices of a list of length val.
    """
    if not isinstance(val, int):
        raise TypeError(_MUST_BE_INT_ERR)
    elif not _is_pow_two_geq_two(val=val):
        raise ValueError(_MUST_BE_POS_INT_ERR)
    elif val not in _REV_IDX_CACHE:
        k: int = val.bit_length() - 1
        _REV_IDX_CACHE[val] = [int(bin(i)[2:].zfill(k)[::-1], 2) for i in range(val)]
    return _REV_IDX_CACHE[val]


def _bit_reverse_copy(val: list) -> list:
    """
    Permute indices by bit-reversal.
    """
    if not isinstance(val, list):
        raise TypeError(_MUST_BE_LIST_ERR)
    n: int = len(val)
    return [val[i] for i in _bit_reverse_len_to_idxs(val=n)]


def _is_rou(root: int, mod: int, deg: int) -> bool:
    """
    Check if val is a root of unity of order 2*deg modulo modulus.
    """
    if (root, mod, deg) not in _ROU_CACHE:
        _ROU_CACHE[(root, mod, deg)] = isinstance(root, int) and isinstance(mod, int) and isinstance(deg, int) and mod >= 3 and deg >= 1 and pow(base=root, exp=2 * deg, mod=mod) == 1
    return _ROU_CACHE[(root, mod, deg)]


def _is_prou(root: int, mod: int, deg: int) -> bool:
    """
    Check if val is a primitive root of order 2*deg modulo modulus.
    """
    if (root, mod, deg) not in _PROU_CACHE:
        is_prim: bool = _is_rou(root=root, mod=mod, deg=deg)
        is_prim = is_prim and all(pow(base=root, exp=i, mod=mod) != 1 for i in range(1, 2*deg))
        _PROU_CACHE[(root, mod, deg)] = is_prim
    return _PROU_CACHE[(root, mod, deg)]


def _find_prou(mod: int, deg: int) -> int:
    """
    Find a primitive root of order 2*deg modulo modulus. Naive loop that first checks 2, then 3, then 4...
    """
    if not isinstance(mod, int) or not isinstance(deg, int):
        raise TypeError(_MUST_BE_INT_ERR)
    elif not _has_prou(mod=mod, deg=deg):
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


def _is_root_inverse(root: int, inv_root: int, mod: int) -> bool:
    return isinstance(root, int) and isinstance(inv_root, int) and isinstance(mod, int) and (root * inv_root) % mod == 1


def _derived_params(deg: int, mod: int, inv_flag: bool) -> _List[int]:
    """ Derive root order from list length, derive root or inv_root from modulus and root order and inv_flag,
    compute root powers, then bit reverse."""
    if (deg, mod, inv_flag) not in _DERIVED_PARAMS_CACHE:
        if not isinstance(deg, int) or not isinstance(mod, int):
            raise TypeError(_MUST_BE_INT_ERR)
        elif not _is_odd_prime(val=mod):
            raise ValueError(_MUST_BE_ODD_PRIME_ERR)
        elif not _has_prou(mod=mod, deg=deg):
            raise ValueError(_MUST_HAVE_PROU_ERR)
        root: int = _find_prou(mod=mod, deg=deg)
        if inv_flag:
            root = pow(base=root, exp=mod-2, mod=mod)
        powers: _List[int] = [pow(base=root, exp=i, mod=mod) for i in range(deg)]
        _DERIVED_PARAMS_CACHE[(deg, mod, inv_flag)] = _bit_reverse_copy(powers)
    return _DERIVED_PARAMS_CACHE[(deg, mod, inv_flag)]


def _is_ntt_valid(val: _List[int], mod: int, inv_flag: bool) -> bool:
    deg: int = len(val)
    brv_powers: _List[int]
    brv_powers = _derived_params(deg=deg, mod=mod, inv_flag=inv_flag)
    val_is_list: bool = isinstance(val, list)
    mod_is_int: bool = isinstance(mod, int)
    inv_flag_is_bool: bool = isinstance(inv_flag, bool)
    brv_powers_is_list: bool = isinstance(brv_powers, list)
    val_and_brv_powers_are_ints: bool = all(isinstance(i, int) for i in val + brv_powers)
    mod_is_odd_prime: bool = _is_odd_prime(val=mod)
    mod_has_prou: bool = _has_prou(mod=mod, deg=deg)
    deg_is_pow_two: bool = _is_pow_two_geq_two(val=deg)
    powers: _List[int] = _bit_reverse_copy(val=brv_powers)
    root: int = powers[1]
    inv_root: int = pow(base=root, exp=mod-2, mod=mod)
    root_is_prou: bool = _is_prou(root=root, mod=mod, deg=deg)
    has_correct_root_inverse: bool = _is_root_inverse(root=root, inv_root=inv_root, mod=mod)
    return all([val_is_list, mod_is_int, inv_flag_is_bool, brv_powers_is_list, val_and_brv_powers_are_ints, mod_is_odd_prime, mod_has_prou, deg_is_pow_two, root_is_prou, has_correct_root_inverse,])


def _cent(val: int, mod: int) -> int:
    """
    Centrally reduce a value modulo q in constant time. Output val satisfies -(q//2) <= val <= q//2.
    """
    if not isinstance(val, int) or not isinstance(mod, int):
        raise TypeError(_MUST_BE_INT_ERR)
    elif mod < 2:
        raise ValueError(_MUST_BE_POS_INT_ERR)
    y: int = val % mod
    w: int = y - (mod//2) - 1
    z: int = y - (1 + (w >> mod.bit_length())) * mod
    return z


def _cooley_tukey_ntt(val: _List[int], mod: int, deg: int, brv_powers: _List[int]) -> _List[int]:
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


def _gentleman_sande_intt(val: _List[int], mod: int, deg: int, brv_powers: _List[int]) -> _List[int]:
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


def _is_ntt_poly_mult_valid(f: _List[int], g: _List[int], deg: int, mod: int, brv_powers: _List[int]):
    f_is_list = isinstance(f, list)
    g_is_list = isinstance(g, list)
    deg_is_int = isinstance(deg, int)
    mod_is_int = isinstance(mod, int)
    brv_powers_is_list = isinstance(brv_powers, list)
    f_and_g_are_ints = all(isinstance(i, int) for i in f+g)
    f_and_g_have_len_deg = len(f) == len(g) == deg
    have_prou = _has_prou(mod=mod, deg=deg)
    powers: _List[int] = _bit_reverse_copy(val=brv_powers)
    root: int = powers[1]
    inv_root: int = pow(base=root, exp=mod - 2, mod=mod)
    have_correct_root = _is_prou(root=root, mod=mod, deg=deg)
    have_correct_inverse_root = _is_root_inverse(root=root, inv_root=inv_root, mod=mod)
    ntt_is_valid_for_f = _is_ntt_valid(val=f, mod=mod, inv_flag=False)
    ntt_is_valid_for_g = _is_ntt_valid(val=g, mod=mod, inv_flag=False)
    return all([f_is_list, g_is_list, deg_is_int, mod_is_int, brv_powers_is_list, f_and_g_are_ints, f_and_g_have_len_deg, have_prou, have_correct_root, have_correct_inverse_root, ntt_is_valid_for_f, ntt_is_valid_for_g])


def _ntt_poly_mult(f: _List[int], g: _List[int], deg: int, mod: int, root: int, inv_root: int, brv_powers: _List[int], brv_inv_root_powers: _List[int]) -> _List[int]:
    """
    Multiply two coefficient representations of polynomials by first computing their NTTs
    and then their Hadamard product before inverting the NTT.
    """
    if not _is_ntt_poly_mult_valid(f=f, g=g, deg=deg, mod=mod, brv_powers=brv_powers):
        raise TypeError(_INVALID_NTT_INPUT_ERR)
    _cooley_tukey_ntt(val=f, mod=mod, deg=deg, brv_powers=brv_powers)
    _cooley_tukey_ntt(val=g, mod=mod, deg=deg, brv_powers=brv_powers)
    fg: _List[int] = [_cent(val=x * y, mod=mod) for x, y in zip(f, g)]
    _gentleman_sande_intt(val=fg, mod=mod, deg=deg, brv_powers=brv_inv_root_powers)
    _gentleman_sande_intt(val=f, mod=mod, deg=deg, brv_powers=brv_inv_root_powers)
    _gentleman_sande_intt(val=g, mod=mod, deg=deg, brv_powers=brv_inv_root_powers)
    return fg
