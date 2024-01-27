from typing import List as _List, Tuple as _Tuple, Dict as _Dict
from copy import deepcopy as _deepcopy
from algebra.ntt import _has_prou, _is_prou, _is_root_inverse, _cent, _derived_params, _cooley_tukey_ntt, _gentleman_sande_intt


def has_prou(mod: int, deg: int) -> bool:
    """ Check if modulus has a primitive root of unity of degree degree"""
    return _has_prou(mod=mod, deg=deg)


def is_prou(root: int, mod: int, deg: int) -> bool:
    """ Check if root is a primitive root of unity mod modulus of degree degree"""
    return _is_prou(root=root, mod=mod, deg=deg)


def is_root_inverse(root: int, inv_root: int, mod: int) -> bool:
    """ Check if root is the inverse of the primitive root of unity mod modulus"""
    return _is_root_inverse(root=root, inv_root=inv_root, mod=mod)


def cent(val: int, mod: int) -> int:
    """ Compute the centered value of val mod modulus"""
    return _cent(val=val, mod=mod)


def ntt(val: _List[int], mod: int, inv_flag: bool) -> _List[int]:
    """ Input a list of integers, a modulus, and a flag indicating forward or inverse, output the NTT"""
    deg: int = len(val)
    brv_powers = _derived_params(deg=deg, mod=mod, inv_flag=inv_flag)
    dcval = _deepcopy(val)
    if not inv_flag:
        return _cooley_tukey_ntt(val=dcval, mod=mod, deg=deg, brv_powers=brv_powers)
    return _gentleman_sande_intt(val=dcval, mod=mod, deg=deg, brv_powers=brv_powers)
