from typing import List as _List, Tuple as _Tuple
from algebra.ntt import _has_prou, _find_prou, _bit_reverse_copy, _cooley_tukey_ntt, _gentleman_sande_intt, _is_odd_prime


def _derived_params(val: _List[int], modulus: int, inv_flag: bool) -> _Tuple[int, int, int, int, int, _List[int]]:
    """ Derive root order from list length, derive root or inv_root from modulus and root order and inv_flag,
    compute root powers, then bit reverse."""
    if not isinstance(val, list) or not all(isinstance(x, int) for x in val) or not isinstance(modulus, int):
        raise TypeError(f"Input must be lists of integers and integer.")
    elif not _is_odd_prime(val=modulus):
        raise ValueError(f"Modulus must be an odd prime")
    elif not _has_prou(mod=modulus, deg=len(val)):
        raise ValueError(f"Modulus does not have primitive root of unity of correct order")
    deg: int = len(val)
    root_order: int = 2*deg
    halfmod: int = modulus//2
    logmod: int = modulus.bit_length()
    root: int = _find_prou(mod=modulus, deg=deg)
    if inv_flag:
        root = pow(base=root, exp=modulus - 2, mod=modulus)
    powers: _List[int] = [pow(base=root, exp=i, mod=modulus) for i in range(deg)]
    brv_powers: _List[int] = _bit_reverse_copy(powers)
    return halfmod, logmod, deg, root_order, root, brv_powers,


def _fwd_ntt_api(val: _List[int], modulus: int) -> _List[int]:
    """ Compute the forward NTT with cooley_tukey_ntt"""
    halfmod, logmod, degree, root_order, root, brv_root_powers = _derived_params(val=val, modulus=modulus, inv_flag=False)
    return _cooley_tukey_ntt(val=val, mod=modulus, halfmod=halfmod, logmod=logmod, deg=root_order, root=root,
                             brv_powers=brv_root_powers)


def _inv_ntt_api(val: _List[int], modulus: int) -> _List[int]:
    """ Compute the inverse NTT with gentleman_sande_intt"""
    halfmod, logmod, degree, root_order, inv_root, brv_inv_root_powers = _derived_params(val=val, modulus=modulus, inv_flag=True)
    return _gentleman_sande_intt(val=val, mod=modulus, halfmod=halfmod, logmod=logmod, deg=root_order, inv_root=inv_root,
                                 brv_powers=brv_inv_root_powers)


def ntt_api(val: _List[int], modulus: int, inv_flag: bool) -> _List[int]:
    """ Input a list of integers, a modulus, and a flag indicating forward or inverse, output the NTT"""
    halfmod, logmod, deg, root_order, root_or_inv_root, powers = _derived_params(val=val, modulus=modulus, inv_flag=inv_flag)
    if not inv_flag:
        return _cooley_tukey_ntt(val=val, mod=modulus, halfmod=halfmod, logmod=logmod, deg=deg, root=root_or_inv_root,
                                 brv_powers=powers)
    return _gentleman_sande_intt(val=val, mod=modulus, halfmod=halfmod, logmod=logmod, deg=deg,
                                 inv_root=root_or_inv_root, brv_powers=powers)
