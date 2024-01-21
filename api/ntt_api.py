from typing import List
from copy import deepcopy
from algebra.ntt import cooley_tukey_ntt, gentleman_sande_intt, derived_params


def fwd_ntt_api(val: List[int], modulus: int) -> List[int]:
    """ Compute the forward NTT with cooley_tukey_ntt"""
    modulus, halfmod, logmod, degree, root_order, root_or_inv_root, brv_powers, inv_flag = derived_params(modulus=modulus, degree=len(val),
                                                                        inv_flag=False)
    return cooley_tukey_ntt(val=deepcopy(val), modulus=modulus, halfmod=halfmod, logmod=logmod, root_order=root_order,
                            bit_rev_root_powers=brv_powers)


def inv_ntt_api(val: List[int], modulus: int) -> List[int]:
    """ Compute the inverse NTT with gentleman_sande_intt"""
    modulus, halfmod, logmod, degree, root_order, root_or_inv_root, brv_powers, inv_flag = derived_params(modulus=modulus, degree=len(val),
                                                                                inv_flag=True)
    return gentleman_sande_intt(val=deepcopy(val), modulus=modulus, halfmod=halfmod, logmod=logmod,
                                root_order=root_order, bit_rev_inv_root_powers=brv_powers)


def ntt_api(val: List[int], modulus: int, inv_flag: bool) -> List[int]:
    """ Input a list of integers, a modulus, and a flag indicating forward or inverse, output the NTT"""
    if not inv_flag:
        return fwd_ntt_api(val=val, modulus=modulus)
    return inv_ntt_api(val=val, modulus=modulus)
