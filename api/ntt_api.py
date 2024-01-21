from typing import List, Tuple
from copy import deepcopy
from algebra.ntt import has_primitive_root_of_unity, find_primitive_root, bit_reverse_copy, cooley_tukey_ntt, gentleman_sande_intt, is_odd_prime


def _derived_params(val: List[int], modulus: int, inv_flag: bool) -> Tuple[int, int, int, int, int, List[int]]:
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


def _fwd_ntt_api(val: List[int], modulus: int) -> List[int]:
    """ Compute the forward NTT with cooley_tukey_ntt"""
    halfmod, logmod, degree, root_order, root, brv_root_powers = _derived_params(val=val, modulus=modulus, inv_flag=False)
    return cooley_tukey_ntt(val=deepcopy(val), modulus=modulus, root_order=root_order, bit_rev_root_powers=brv_root_powers, halfmod=halfmod, logmod=logmod)


def _inv_ntt_api(val: List[int], modulus: int) -> List[int]:
    """ Compute the inverse NTT with gentleman_sande_intt"""
    halfmod, logmod, degree, root_order, inv_root, brv_inv_root_powers = _derived_params(val=val, modulus=modulus, inv_flag=True)
    return gentleman_sande_intt(val=deepcopy(val), modulus=modulus, root_order=root_order, bit_rev_inv_root_powers=brv_inv_root_powers, halfmod=halfmod, logmod=logmod)


def ntt_api(val: List[int], modulus: int, inv_flag: bool) -> List[int]:
    """ Input a list of integers, a modulus, and a flag indicating forward or inverse, output the NTT"""
    if not inv_flag:
        return _fwd_ntt_api(val=val, modulus=modulus)
    return _inv_ntt_api(val=val, modulus=modulus)
