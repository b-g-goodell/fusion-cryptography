from typing import List
from copy import deepcopy
from algebra.ntt import cooley_tukey_ntt, gentleman_sande_intt, derived_params


def fwd_ntt_api(val: List[int], modulus: int) -> List[int]:
    """ Compute the forward NTT with cooley_tukey_ntt"""
    halfmod, logmod, degree, root_order, root, brv_root_powers = derived_params(degree=len(val), modulus=modulus, inv_flag=False)
    return cooley_tukey_ntt(val=deepcopy(val), modulus=modulus, root_order=root_order, bit_rev_root_powers=brv_root_powers, halfmod=halfmod, logmod=logmod)


def inv_ntt_api(val: List[int], modulus: int) -> List[int]:
    """ Compute the inverse NTT with gentleman_sande_intt"""
    halfmod, logmod, degree, root_order, inv_root, brv_inv_root_powers = derived_params(degree=len(val), modulus=modulus, inv_flag=True)
    return gentleman_sande_intt(val=deepcopy(val), modulus=modulus, root_order=root_order, bit_rev_inv_root_powers=brv_inv_root_powers, halfmod=halfmod, logmod=logmod)


def ntt_api(val: List[int], modulus: int, inv_flag: bool) -> List[int]:
    """ Input a list of integers, a modulus, and a flag indicating forward or inverse, output the NTT"""
    if not inv_flag:
        return fwd_ntt_api(val=val, modulus=modulus)
    return inv_ntt_api(val=val, modulus=modulus)
