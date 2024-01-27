from typing import List as _List
from copy import deepcopy as _deepcopy
from algebra.ntt import (_find_prou, _bit_reverse_copy, _has_prou, _is_prou,
                         _is_root_inverse, _cent, _brv_root_and_inv_root_powers, _cooley_tukey_ntt,
                         _gentleman_sande_intt, _ntt_poly_mult)


def bit_reverse_copy(val: list) -> list:
    """ Input a list, output a bit reversed copy of the list"""
    return _bit_reverse_copy(val=val)


def find_prou(mod: int, deg: int) -> int:
    """ Find a primitive root of unity of degree 2*degree mod modulus"""
    return _find_prou(mod=mod, deg=deg)


def has_prou(mod: int, deg: int) -> bool:
    """ Check if modulus has a primitive root of unity of degree 2*degree"""
    return _has_prou(mod=mod, deg=deg)


def is_prou(root: int, mod: int, deg: int) -> bool:
    """ Check if root is a primitive root of unity mod modulus of degree 2*degree"""
    return _is_prou(root=root, mod=mod, deg=deg)


def is_root_inverse(root: int, inv_root: int, mod: int) -> bool:
    """ Check if root is the inverse of the primitive root of unity mod modulus"""
    return _is_root_inverse(root=root, inv_root=inv_root, mod=mod)


def cent(val: int, mod: int) -> int:
    """ Compute the centered value of val mod modulus"""
    return _cent(val=val, mod=mod)


def ntt_poly_mult(f: _List[int], g: _List[int], mod: int) -> _List[int]:
    """ Input two lists of integers and a modulus, output the NTT of their product"""
    deg: int = len(f)
    brv_powers: _List[int]
    brv_inv_powers: _List[int]
    brv_powers, brv_inv_powers = _brv_root_and_inv_root_powers(deg=deg, mod=mod)
    root: int
    inv_root: int
    root, inv_root = bit_reverse_copy(val=brv_powers)[1], bit_reverse_copy(val=brv_inv_powers)[1]
    dcf: _List[int]
    dcg: _List[int]
    dcf, dcg = _deepcopy(f), _deepcopy(g)
    return _ntt_poly_mult(f=dcf, g=dcg, mod=mod, brv_powers=brv_powers, brv_inv_root_powers=brv_inv_powers)


def ntt(val: _List[int], mod: int, inv_flag: bool) -> _List[int]:
    """ Input a list of integers, a modulus, and a flag indicating forward or inverse, output the NTT"""
    deg: int = len(val)
    brv_powers: _List[int]
    brv_inv_powers: _List[int]
    brv_powers, brv_inv_powers = _brv_root_and_inv_root_powers(deg=deg, mod=mod)
    dcval = _deepcopy(val)
    if not inv_flag:
        return _cooley_tukey_ntt(val=dcval, mod=mod, brv_powers=brv_powers)
    return _gentleman_sande_intt(val=dcval, mod=mod, brv_powers=brv_inv_powers)
