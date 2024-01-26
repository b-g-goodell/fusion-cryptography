from typing import List as _List, Tuple as _Tuple
from copy import deepcopy
from algebra.ntt import _has_prou, _find_prou, _bit_reverse_copy, _cooley_tukey_ntt, _gentleman_sande_intt, _is_odd_prime, _cent, _is_prou, _is_root_inverse
from algebra.errors import _MUST_BE_LIST_ERR, _MUST_BE_ODD_PRIME_ERR, _MUST_HAVE_PROU_ERR


def _derived_params(val: _List[int], mod: int, inv_flag: bool) -> _Tuple[int, int, int, int, int, _List[int]]:
    """ Derive root order from list length, derive root or inv_root from modulus and root order and inv_flag,
    compute root powers, then bit reverse."""
    if not isinstance(val, list):
        raise TypeError(_MUST_BE_LIST_ERR)
    elif not _is_odd_prime(val=mod):
        raise ValueError(_MUST_BE_ODD_PRIME_ERR)
    elif not _has_prou(mod=mod, deg=len(val)):
        raise ValueError(_MUST_HAVE_PROU_ERR)
    deg: int = len(val)
    root_order: int = 2*deg
    halfmod: int = mod // 2
    logmod: int = mod.bit_length()
    root: int = _find_prou(mod=mod, deg=deg)
    if inv_flag:
        root = pow(base=root, exp=mod - 2, mod=mod)
    powers: _List[int] = [pow(base=root, exp=i, mod=mod) for i in range(deg)]
    brv_powers: _List[int] = _bit_reverse_copy(powers)
    return halfmod, logmod, deg, root_order, root, brv_powers,


def cent(val: int, mod: int) -> int:
    """ Compute the centered value of val mod modulus"""
    return _cent(val=val, mod=mod)


def is_prou(root: int, mod: int, deg: int) -> bool:
    """ Check if root is a primitive root of unity mod modulus of degree degree"""
    return _is_prou(root=root, mod=mod, deg=deg)


def has_prou(mod: int, deg: int) -> bool:
    """ Check if modulus has a primitive root of unity of degree degree"""
    return _has_prou(mod=mod, deg=deg)


def is_root_inverse(root: int, inv_root: int, mod: int) -> bool:
    """ Check if root is the inverse of the primitive root of unity mod modulus"""
    return _is_root_inverse(root=root, inv_root=inv_root, mod=mod)


def _fwd_ntt(val: _List[int], mod: int) -> _List[int]:
    """ Compute the forward NTT with cooley_tukey_ntt"""
    halfmod, logmod, degree, root_order, root, brv_root_powers = _derived_params(val=val, mod=mod, inv_flag=False)
    return _cooley_tukey_ntt(val=val, mod=mod, deg=root_order, root=root, brv_powers=brv_root_powers)


def _inv_ntt(val: _List[int], mod: int) -> _List[int]:
    """ Compute the inverse NTT with gentleman_sande_intt"""
    halfmod, logmod, degree, root_order, inv_root, brv_inv_root_powers = _derived_params(val=val, mod=mod,
                                                                                         inv_flag=True)
    return _gentleman_sande_intt(val=val, mod=mod, deg=root_order, inv_root=inv_root,
                                 brv_powers=brv_inv_root_powers)


def ntt(val: _List[int], mod: int, inv_flag: bool) -> _List[int]:
    """ Input a list of integers, a modulus, and a flag indicating forward or inverse, output the NTT"""
    halfmod, logmod, deg, root_order, root_or_inv_root, powers = _derived_params(val=val, mod=mod,
                                                                                 inv_flag=inv_flag)
    dcval = deepcopy(val)
    if not inv_flag:
        return _cooley_tukey_ntt(val=dcval, mod=mod, deg=deg, root=root_or_inv_root, brv_powers=powers)
    return _gentleman_sande_intt(val=dcval, mod=mod, deg=deg, inv_root=root_or_inv_root, brv_powers=powers)
