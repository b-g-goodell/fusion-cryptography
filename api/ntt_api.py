# Let's mock up a simple API
from typing import List, Tuple
from algebra.ntt import has_primitive_root_of_unity, find_primitive_root, bit_reverse_copy, cooley_tukey_ntt, gentleman_sande_intt, is_odd_prime


def _val_and_modulus_to_halfmod_logmod_degree_and_root_order(val: List[int], modulus: int):
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
    return halfmod, logmod, degree, root_order


def _expand_val_and_modulus_for_fwd_ntt(val: List[int], modulus: int) -> Tuple[List[int], int, int, int, int, int, int, List[int]]:
    halfmod, logmod, degree, root_order = _val_and_modulus_to_halfmod_logmod_degree_and_root_order(val=val, modulus=modulus)
    root: int = find_primitive_root(modulus=modulus, root_order=root_order)
    root_powers: List[int] = [pow(base=root, exp=i, mod=modulus) for i in range(root_order)]
    brv_root_powers: List[int] = bit_reverse_copy(root_powers)
    return halfmod, logmod, degree, root_order, root, brv_root_powers,


def _expand_val_and_modulus_for_inv_ntt(val: List[int], modulus: int) -> Tuple[List[int], int, int, int, int, int, int, List[int]]:
    halfmod, logmod, degree, root_order = _val_and_modulus_to_halfmod_logmod_degree_and_root_order(val=val, modulus=modulus)
    root: int = find_primitive_root(modulus=modulus, root_order=root_order)
    inv_root: int = pow(base=root, exp=modulus - 2, mod=modulus)
    inv_root_powers: List[int] = [pow(base=inv_root, exp=i, mod=modulus) for i in range(root_order)]
    brv_inv_root_powers: List[int] = bit_reverse_copy(inv_root_powers)
    return halfmod, logmod, degree, root_order, inv_root, brv_inv_root_powers,


def _fwd_ntt_api(val: List[int], modulus: int) -> List[int]:
    halfmod, logmod, degree, root_order, root, brv_root_powers = _expand_val_and_modulus_for_fwd_ntt(val=val, modulus=modulus)
    return cooley_tukey_ntt(val=val, modulus=modulus, root_order=root_order, bit_rev_root_powers=brv_root_powers, halfmod=halfmod, logmod=logmod)


def _inv_ntt_api(val: List[int], modulus: int) -> List[int]:
    halfmod, logmod, degree, root_order, inv_root, brv_inv_root_powers = _expand_val_and_modulus_for_inv_ntt(val=val, modulus=modulus)
    return gentleman_sande_intt(val=val, modulus=modulus, root_order=root_order, bit_rev_inv_root_powers=brv_inv_root_powers, halfmod=halfmod, logmod=logmod)


def ntt_api(val: List[int], modulus: int, inv_flag: bool) -> List[int]:
    if not inv_flag:
        return _fwd_ntt_api(val=val, modulus=modulus)
    return _inv_ntt_api(val=val, modulus=modulus)
