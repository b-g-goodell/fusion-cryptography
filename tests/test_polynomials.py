import pytest
from copy import deepcopy
from random import randrange
from typing import List, Tuple
from algebra.polynomials import (
    _validate_modulus_and_vals,
    PolynomialCoefficientRepresentation as Poly,
    PolynomialNTTRepresentation as PolyNTT,
    transform,
)
from algebra.ntt import cent, bit_reverse_copy, ntt_poly_mult, cooley_tukey_ntt, gentleman_sande_intt, find_primitive_root, derived_params
from test_ntt import SAMPLE_SIZE, PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE


RING_PARAMETERS: List[tuple] = [
    # d q q q//2, log(q) d 2d, s, [s], flag
    t + derived_params(modulus=t[1], degree=t[0], inv_flag=False)
    for t in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE
]
# q q//2, log(q) d 2d, s, [s], flag
RING_PARAMETERS: List[tuple] = [t[2:] for t in RING_PARAMETERS]

VALIDATE_MODULUS_AND_COEFFICIENTS_TEST_DATA = [
    t + tuple([
        list(range(t[3]))  # dummy list of correct length
    ])
    for t in RING_PARAMETERS
]
print()

@pytest.mark.parametrize("modulus,halfmod,logmod,degree,root_order,root,brv_powers,inv_flag,dummy_val", VALIDATE_MODULUS_AND_COEFFICIENTS_TEST_DATA)
def test_validate_modulus_and_coefficients(modulus,halfmod,logmod,degree,root_order,root,brv_powers,inv_flag,dummy_val):
    _validate_modulus_and_vals(modulus=modulus, vals=dummy_val)


VALIDATE_MODULUS_AND_COEFFICIENTS_ERRORS = [
    ("not an int", "not a list of ints", TypeError, "modulus must be an int"),
    (17, "not a list of ints", TypeError, "coefficients must be list of integers"),
    (17, ["not a list of ints"] * 3, TypeError, "coefficients must be list of integers"),
    (16, list(range(3)), ValueError, "modulus must be an odd prime"),
    (17, list(range(3)), ValueError, "modulus does not have primitive root of unity of appropriate order."),
    (97, list(range(48)), ValueError, "coefficient vector must have power of 2 length"),
]
print()

@pytest.mark.parametrize("modulus,coefs,exception_type,error_message", VALIDATE_MODULUS_AND_COEFFICIENTS_ERRORS)
def test_validate_modulus_and_coefficients_errors(modulus,coefs,exception_type,error_message):
    with pytest.raises(exception_type) as exc_info:
        _validate_modulus_and_vals(modulus=modulus, vals=coefs)
    assert str(exc_info.value) == error_message


VALIDATE_SAME_RING_TEST_DATA = [t + s for t in VALIDATE_MODULUS_AND_COEFFICIENTS_TEST_DATA for s in VALIDATE_MODULUS_AND_COEFFICIENTS_TEST_DATA]


@pytest.mark.parametrize("degree,modulus,halfmod,logmod,root_order,root,brv_powers,coefs", VALIDATE_SAME_RING_TEST_DATA)
def test_validate_same_ring(degree,modulus,halfmod,logmod,root_order,root,brv_powers,coefs):
    _validate_modulus_and_vals(modulus=modulus, vals=coefs)
#
#
# VALIDATE_MODULUS_AND_COEFFICIENTS_ERRORS = [
#     ("not an int", "not a list of ints", TypeError, "modulus must be an int"),
#     (17, "not a list of ints", TypeError, "coefficients must be list of integers"),
#     (17, ["not a list of ints"] * 3, TypeError, "coefficients must be list of integers"),
#     (16, list(range(3)), ValueError, "modulus must be an odd prime"),
#     (17, list(range(3)), ValueError, "modulus does not have primitive root of unity of appropriate order."),
#     (97, list(range(48)), ValueError, "coefficient vector must have power of 2 length"),
# ]
# print()
#
# @pytest.mark.parametrize("modulus,coefs,exception_type,error_message", VALIDATE_MODULUS_AND_COEFFICIENTS_ERRORS)
# def test_validate_modulus_and_coefficients_errors(modulus,coefs,exception_type,error_message):
#     with pytest.raises(exception_type) as exc_info:
#         _validate_modulus_and_vals(modulus=modulus, vals=coefs)
#     assert str(exc_info.value) == error_message
#
#
# # CHECK_MODULUS_HALFMOD_LOGMOD_TEST_DATA = [
# #     (x, x // 2, x.bit_length())  # adjusted logmod to be x.bit_length() - 1
# #     for x in range(3, 100)
# # ]
# # CHECK_MODULUS_HALFMOD_LOGMOD_ERRORS = [
# #     ("not an int", 3, 3, TypeError, "Modulus, halfmod, and logmod must all be integers."),
# #     (6, "not an int", 3, TypeError, "Modulus, halfmod, and logmod must all be integers."),
# #     (6, 3, "not an int", TypeError, "Modulus, halfmod, and logmod must all be integers."),
# #     (2, 1, 2, ValueError, "Modulus must be >=3, had 2"),
# #     (3, 0, 2, ValueError, "Halfmod must be >=1, had 0"),
# #     (3, 2, 2, ValueError, "Halfmod must be half modulus, but had halfmod=2 and modulus=3"),
# #     (6, 2, 0, ValueError, "Halfmod must be half modulus, but had halfmod=2 and modulus=6"),
# #     (6, 3, 2, ValueError,
# #      "Logmod must bound the modulus, had 2**(logmod-1)=2, 2**logmod=4, and modulus=6")  # adjust the error message to match the actual output
# #     # The above line has been wrapped for readability, but should be a single line in the code.
# # ]
#
# # @pytest.mark.parametrize("modulus,halfmod,logmod", CHECK_MODULUS_HALFMOD_LOGMOD_TEST_DATA)
# # def test_check_modulus_halfmod_logmod_valid(modulus, halfmod, logmod):
# #     check_modulus_halfmod_logmod(modulus, halfmod, logmod)
# #
# #
# # @pytest.mark.parametrize("modulus,halfmod,logmod,exception_type,error_message", CHECK_MODULUS_HALFMOD_LOGMOD_ERRORS)
# # def test_check_modulus_halfmod_logmod_invalid(modulus, halfmod, logmod, exception_type, error_message):
# #     with pytest.raises(exception_type) as exc_info:
# #         check_modulus_halfmod_logmod(modulus, halfmod, logmod)
# #     assert str(exc_info.value) == error_message
