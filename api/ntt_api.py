from typing import List as _List
from copy import deepcopy as _deepcopy
from algebra.ntt import (
    _cooley_tukey_ntt, _gentleman_sande_intt, _derived_params, _cent, _find_primitive_root,
    _has_primitive_root_of_unity, _is_odd_prime, _is_pow_two_geq_two, _ERR_MUST_BE_INT,
    _ERR_MUST_BE_LIST, _ERR_MUST_BE_BOOL, _ERR_MUST_BE_POW_TWO, _ERR_MUST_BE_GEQ_THREE, _ERR_MUST_BE_STRICTLY_POSITIVE,
    _ERR_MUST_BE_FLOOR_HALF, _ERR_MUST_BE_BIT_LENGTH, _ERR_MUST_HAVE_PRIM_ROU, _ERR_MUST_BE_ODD_PRIME,
    _ERR_MUST_BE_CORRECT_ROOT, _ERR_MUST_BE_SAME_LENGTH,
)

ERR_MUST_BE_INT: str = _ERR_MUST_BE_INT
ERR_MUST_BE_LIST: str = _ERR_MUST_BE_LIST
ERR_MUST_BE_BOOL: str = _ERR_MUST_BE_BOOL
ERR_MUST_BE_POW_TWO: str = _ERR_MUST_BE_POW_TWO
ERR_MUST_BE_GEQ_THREE: str = _ERR_MUST_BE_GEQ_THREE
ERR_MUST_BE_STRICTLY_POSITIVE: str = _ERR_MUST_BE_STRICTLY_POSITIVE
ERR_MUST_BE_FLOOR_HALF: str = _ERR_MUST_BE_FLOOR_HALF
ERR_MUST_BE_BIT_LENGTH: str = _ERR_MUST_BE_BIT_LENGTH
ERR_MUST_HAVE_PRIM_ROU: str = _ERR_MUST_HAVE_PRIM_ROU
ERR_MUST_BE_ODD_PRIME: str = _ERR_MUST_BE_ODD_PRIME
ERR_MUST_BE_CORRECT_ROOT: str = _ERR_MUST_BE_CORRECT_ROOT
ERR_MUST_BE_SAME_LENGTH: str = _ERR_MUST_BE_SAME_LENGTH


def cent(val: int, mod: int) -> int:
    return _cent(val=val, modulus=mod, halfmod=mod//2, logmod=mod.bit_length())


def find_primitive_root(mod: int, root_order: int) -> int:
    return _find_primitive_root(modulus=mod, root_order=root_order)


def has_primitive_root_of_unity(mod: int, root_order: int) -> bool:
    return _has_primitive_root_of_unity(modulus=mod, root_order=root)


def is_odd_prime(val: int) -> bool:
    return _is_odd_prime(val=val)


def is_pow_two_geq_two(val: int) -> bool:
    return _is_pow_two_geq_two(val=val)


def fwd_ntt_api(val: _List[int], modulus: int) -> _List[int]:
    """ Compute the forward NTT with cooley_tukey_ntt"""
    modulus, halfmod, logmod, degree, root_order, root_or_inv_root, brv_powers, inv_flag = _derived_params(modulus=modulus, degree=len(val),
                                                                                                           inv_flag=False)
    return _cooley_tukey_ntt(val=_deepcopy(val), modulus=modulus, halfmod=halfmod, logmod=logmod, root_order=root_order,
                             bit_rev_root_powers=brv_powers)


def inv_ntt_api(val: _List[int], modulus: int) -> _List[int]:
    """ Compute the inverse NTT with gentleman_sande_intt"""
    modulus, halfmod, logmod, degree, root_order, root_or_inv_root, brv_powers, inv_flag = _derived_params(modulus=modulus, degree=len(val),
                                                                                                           inv_flag=True)
    return _gentleman_sande_intt(val=_deepcopy(val), modulus=modulus, halfmod=halfmod, logmod=logmod,
                                 root_order=root_order, bit_rev_inv_root_powers=brv_powers)


def ntt_api(val: _List[int], modulus: int, inv_flag: bool) -> _List[int]:
    """ Input a list of integers, a modulus, and a flag indicating forward or inverse, output the NTT"""
    if not inv_flag:
        return fwd_ntt_api(val=val, modulus=modulus)
    return inv_ntt_api(val=val, modulus=modulus)
