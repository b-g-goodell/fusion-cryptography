# import pytest
# from math import ceil, log2
# from typing import List, Dict, Union
# from algebra.sampling import sample_coefficient_matrix, sample_ntt_matrix
# from algebra.polynomials import (
#     PolynomialCoefficientRepresentation as Poly,
#     PolynomialNTTRepresentation as PolyNTT,
#     transform,
# )
# from fusion.fusion import (
#     PREFIX_PARAMETERS,
#     fusion_setup,
#     fusion_keygen,
#     fusion_sign,
#     # aggregate,
#     # verify,
#     Params,
#     OneTimeKeyTuple,
#     OneTimeVerificationKey,
#     Signature,
#     OneTimeSigningKey,
#     GeneralMatrix,
#     SignatureChallenge,
#     _pre_hash_msg_to_int_digest,
#     _hash_vk_and_pre_hashed_msg_to_bytes_digest,
#     AggregationCoefficient,
#     # _make_agg_coefs,
#     _bytes_digest_to_polys,
#     sha3_256_wrapper,
#     shake_256_wrapper,
#     _make_sig_challenge,
#     bits_for_bdd_coef,
#     bits_for_index,
#     bits_for_fy_shuffle,
#     bits_for_bdd_poly,
#     bytes_for_bdd_coef,
#     bytes_for_bdd_poly,
#     bytes_for_fy_shuffle,
#     bytes_for_index
# )
#
# TEST_SAMPLE_SIZE: int = 2**5
#
# # Test global constants
#
# def test_prefix_parameters():
#     expected_param_keys = [
#         "capacity", "modulus", "degree", "root_order", "root", "inv_root",
#         "num_rows_pub_challenge", "num_rows_sk", "num_rows_vk", "num_cols_sk",
#         "num_cols_vk", "sign_pre_hash_dst", "sign_hash_dst", "agg_xof_dst",
#         "beta_sk", "beta_ch", "beta_ag", "omega_sk", "omega_ch", "omega_ag",
#         "omega_vf_intermediate", "beta_vf_intermediate", "omega_vf", "beta_vf",
#         # "bytes_for_one_coef_bdd_by_beta_ch", "bytes_for_one_coef_bdd_by_beta_ag",
#         # "bytes_for_index", "bytes_for_fy_shuffle", "bytes_per_sig_chall_poly",
#         # "bytes_per_agg_coef_poly"
#     ]
#     assert 128 in PREFIX_PARAMETERS
#     assert 256 in PREFIX_PARAMETERS
#     for secpar in [128, 256]:
#         these_prefixed_parameters = PREFIX_PARAMETERS[secpar]
#         assert all(next_param_key in these_prefixed_parameters for next_param_key in expected_param_keys)
#
# # Test simple wrappers for sha3 and shake
#
# def dummy_sha3_256(input_data: bytes) -> bytes:
#     return b'\x00\x01' * 16
#
#
# def dummy_shake_256(input_data: bytes, n: int):
#     # Return a fixed byte sequence of length n for testing
#     return b'\xab\xcd' * n
#
#
# def test_sha3_256_wrapper():
#     # The expected hash output from the sha3_256 algorithm for "Hello, World!" in bytes
#     expected_hash = bytes.fromhex("1af17a664e3fa8e419b8ba05c2a173169df76162a5a286e0c405b460d478f7ef")
#
#     # The result of the sha3_256_wrapper function
#     result = sha3_256_wrapper("Hello, World!")
#
#     # Assert that the expected hash matches the calculated hash
#     assert result == expected_hash, f"Expected {expected_hash.hex()}, but got {result.hex()}"
#
#
# def test_shake_256_wrapper():
#     # The expected output from the shake_256 algorithm for "Hello, World!", requesting 512 bits (64 bytes)
#     expected_output = bytes.fromhex(
#         "b3be97bfd978833a65588ceae8a34cf59e95585af62063e6b89d0789f372424e8b0d1be4f21b40ce5a83a438473271e0661854f02d431db74e6904d6c347d757")
#     # The result of the shake_256_wrapper function
#     result = shake_256_wrapper("Hello, World!", 64)
#     # Assert that the expected output matches the calculated output
#     assert result == expected_output, f"Expected {expected_output.hex()}, but got {result.hex()}"
#
# # Test simple functions for counting bits and bytes
#
# @pytest.mark.parametrize(
#     "secpar, beta, expected_bits, expected_bytes",
#     [
#         (128, 1, 130, 17),    # Simple case
#         (0, 0, 0, 0),          # Boundary case where secpar and beta are 0
#         (128, 0, 128, 16),     # Case with beta as 0
#         (0, 1, 2, 1),          # Case with secpar as 0
#         (256, 1000, 267, 34),  # Larger secpar and beta values
#     ]
# )
# def test_bits_and_bytes_for_bdd_coef(secpar, beta, expected_bits, expected_bytes):
#     # Test bits calculation
#     calculated_bits = bits_for_bdd_coef(secpar=secpar, beta=beta)
#     calculated_bytes = bytes_for_bdd_coef(secpar=secpar, beta=beta)
#     assert calculated_bits == expected_bits, f"Expected bits for secpar {secpar} and beta {beta} to be {expected_bits}, got {calculated_bits}"
#     assert calculated_bytes == expected_bytes, f"Expected bytes for secpar {secpar} and beta {beta} to be {expected_bytes}, got {calculated_bytes}"
#     assert ceil(expected_bits/8) == expected_bytes, f"Expected bytes for secpar {secpar} and beta {beta} to be {expected_bytes}, got {ceil(expected_bits/8)}"
#
#
# @pytest.mark.parametrize(
#     "secpar, degree, expected_bits, expected_bytes",
#     [
#         (128, 2, 129, 17),    # Simple case
#         (0, 2, 1, 1),          # Case with secpar as 0
#         (256, 1000, 266, 34),  # Larger secpar and beta values
#         # Add more cases including failure
#     ]
# )
# def test_bits_and_bytes_for_index(secpar, degree, expected_bits, expected_bytes):
#     # Test bits calculation
#     calculated_bits = bits_for_index(secpar=secpar, degree=degree)
#     calculated_bytes = bytes_for_index(secpar=secpar, degree=degree)
#     assert calculated_bits == expected_bits, f"Expected bits for secpar {secpar} and degree {degree} to be {expected_bits}, got {calculated_bits}"
#     assert calculated_bytes == expected_bytes, f"Expected bytes for secpar {secpar} and degree {degree} to be {expected_bytes}, got {calculated_bytes}"
#     assert ceil(expected_bits/8) == expected_bytes, f"Expected bytes for secpar {secpar} and degree {degree} to be {expected_bytes}, got {ceil(expected_bits/8)}"
#
#
# @pytest.mark.parametrize(
#     "secpar, degree, expected_bits, expected_bytes",
#     [
#         (128, 2, 258, 33),    # Simple case
#         (0, 2, 2, 1),          # Case with secpar as 0
#         (256, 1000, 266000, 33250),  # Larger secpar and beta values
#         # Add more cases as necessary
#     ]
# )
# def test_bits_and_bytes_for_fy_shuffle(secpar, degree, expected_bits, expected_bytes):
#     calculated_bits = bits_for_fy_shuffle(secpar=secpar, degree=degree)
#     calculated_bytes = bytes_for_fy_shuffle(secpar=secpar, degree=degree)
#     assert calculated_bits == expected_bits, f"Expected bits for secpar {secpar} and degree {degree} to be { expected_bits}, got { calculated_bits}"
#     assert calculated_bytes == expected_bytes, f"Expected bytes for secpar {secpar} and degree {degree} to be {expected_bytes}, got { calculated_bytes}"
#     assert ceil(expected_bits/8) == expected_bytes
#
#
# @pytest.mark.parametrize(
#     "secpar, beta, degree, omega, expected_bits, expected_bytes",
#     [
#         (128, 1, 8, 6, 1828, 229),    # Simple case
#         (0, 0, 8, 0, 24, 3),  # An edge case
#         (0, 0, 4, 0, 8, 1),  # An edge case
#         (0, 0, 2, 0, 2, 1),  # An edge case
#         (0, 0, 8, 1, 24, 3),  # An edge case
#         (0, 1, 8, 0, 24, 3),  # An edge case
#         (1, 0, 8, 0, 32, 4),  # An edge case
#         (0, 1, 8, 1, 26, 4),  # An edge case
#         (1, 0, 8, 1, 33, 5),  # An edge case
#         (1, 1, 8, 0, 32, 4),  # An edge case
#         (1, 1, 8, 1, 35, 5),  # small non-edge case
#         (256, 1000, 1000, 50, 279350, 34919)
#     ]
# )
# def test_bits_and_bytes_for_bdd_poly(secpar, beta, degree, omega, expected_bits, expected_bytes):
#     calculated_bits = bits_for_bdd_poly( secpar=secpar, beta=beta, degree=degree, omega=omega)
#     calculated_bytes = bytes_for_bdd_poly( secpar=secpar, beta=beta, degree=degree, omega=omega)
#     assert calculated_bits == expected_bits, (f"Expected bits for secpar {secpar} and degree {degree} to be { expected_bits}, got { calculated_bits}")
#     assert calculated_bytes == expected_bytes, (f"Expected bytes for secpar {secpar} and degree {degree} to be { expected_bytes}, got { calculated_bytes}")
#     assert calculated_bytes == ceil(calculated_bits/8), (f"Expected bytes for secpar {secpar} and degree {degree} to be { ceil(calculated_bits)} and got { calculated_bytes}")
#
#
# # Test classes for typing
#
# # Test Params class
# def test_params_initialization_success():
#     # Assuming PREFIX_PARAMETERS contains at least one secpar key
#     test_secpar = list(PREFIX_PARAMETERS.keys())[0]
#     params = Params(test_secpar)
#     attributes = PREFIX_PARAMETERS[test_secpar]
#
#     for key, value in attributes.items():
#         assert getattr(params, key) == value
#
#
# def test_params_initialization_failure():
#     with pytest.raises(ValueError):
#         Params('invalid_secpar')
#
#
# @pytest.mark.parametrize('secpar', list(PREFIX_PARAMETERS.keys()))
# def test_assigned_attributes(secpar):
#     params = Params(secpar)
#     expected_values = PREFIX_PARAMETERS[secpar]
#
#     for attr, expected in expected_values.items():
#         assert getattr(params, attr) == expected
#
#
# @pytest.mark.parametrize('secpar', list(PREFIX_PARAMETERS.keys()))
# def test_str_representation(secpar):
#     params = Params(secpar)
#     assert str(params).startswith('Params(')
#
#
# @pytest.mark.parametrize('secpar', list(PREFIX_PARAMETERS.keys()))
# def test_repr_representation(secpar):
#     params = Params(secpar)
#     assert repr(params) == str(params)
#
#
# def test_params_equality():
#     test_secpar = list(PREFIX_PARAMETERS.keys())[0]
#     params1 = Params(test_secpar)
#     params2 = Params(test_secpar)
#     assert params1 == params2
#     assert params1 is not params2
#
#
# def test_params_inequality():
#     test_secpars = list(PREFIX_PARAMETERS.keys())[:2]
#     assert test_secpars
#     assert len(test_secpars) == 2 # Ensure there are at least two elements
#     params1 = Params(test_secpars[0])
#     params2 = Params(test_secpars[1])
#     assert params1 != params2
#
#
# # Test public_challenge
# @pytest.mark.parametrize('secpar', list(PREFIX_PARAMETERS.keys()))
# def test_params_public_challenge(secpar):
#     params = Params(secpar)
#     assert isinstance(params.public_challenge, GeneralMatrix)
#     assert params.public_challenge.elem_class == PolyNTT
#     assert len(params.public_challenge.matrix) == params.num_rows_pub_challenge
#     assert all(len(y) == params.num_cols_pub_challenge for y in params.public_challenge.matrix)
#     halfmod = params.modulus//2
#     assert all(z.norm() <= halfmod for y in params.public_challenge.matrix for z in y)
#     assert all(z.weight() <= params.degree for y in params.public_challenge.matrix for z in y)
#
#
# # Test otsk class
# # Mock an algebraic class with minimal implementation for testing GeneralMatrix.
# class TestAlgebraicClass:
#     def __init__(self, value):
#         self.value = value
#
#     def __add__(self, other):
#         return TestAlgebraicClass(self.value + other.value)
#
#     def __radd__(self, other):
#         return self + other
#
#     def __neg__(self):
#         return TestAlgebraicClass(-self.value)
#
#     def __sub__(self, other):
#         return TestAlgebraicClass(self.value - other.value)
#
#     def __mul__(self, other):
#         if not isinstance(other, TestAlgebraicClass):
#             return NotImplemented
#         return TestAlgebraicClass(self.value * other.value)
#
#     def __eq__(self, other):
#         if not isinstance(other, TestAlgebraicClass):
#             return NotImplemented
#         return self.value == other.value
#
#     def __mod__(self, other):
#         return TestAlgebraicClass(self.value % other)
#
#     def norm(self, p: Union[int, str] = "infty"):
#         return abs(self.value)
#
#     def weight(self):
#         return self.value
#
#
# # Function to create a GeneralMatrix instance with mock TestAlgebraicClass instances.
# def create_general_matrix(rows, cols):
#     return GeneralMatrix([[TestAlgebraicClass(i + j) for j in range(cols)] for i in range(rows)])
#
#
# @pytest.fixture
# def general_matrix_2x2():
#     # 2x2 GeneralMatrix fixture for testing
#     return create_general_matrix(2, 2)
#
#
# @pytest.fixture
# def general_matrix_3x2():
#     # 3x2 GeneralMatrix fixture for testing (to test differing sizes)
#     return create_general_matrix(3, 2)
#
#
# # Mock the is_algebraic_class function to always return True for the test algebraic class.
#
# @pytest.fixture(autouse=True)
# def mock_is_algebraic_class(monkeypatch):
#     monkeypatch.setattr(
#
#         'algebra.matrices.is_algebraic_class',
#
#         lambda cls: cls == TestAlgebraicClass
#
#     )
#
#
# def is_algebraic_class(cls):
#     return True
#
#
# def test_one_time_signing_key_initialization():
#     matrix1 = GeneralMatrix([[1, 2], [3, 4]])
#     matrix2 = GeneralMatrix([[5, 6], [7, 8]])
#
#     # Test valid initialization
#     signing_key = OneTimeSigningKey(matrix1, matrix2)
#
#     # Test whether string representation is as expected
#     assert str(signing_key) == (
#         "OneTimeSigningKey(left_sk_hat=GeneralMatrix(elem_class=<class 'int'>, matrix=[[1, 2], [3, 4]]), "
#         "right_sk_hat=GeneralMatrix(elem_class=<class 'int'>, matrix=[[5, 6], [7, 8]]))"), "String representation of the key is not correct."
#
#     # Test invalid initialization with mismatched matrix dimensions which should raise ValueError
#     matrix3 = GeneralMatrix([[1, 2, 3], [4, 5, 6]])
#     with pytest.raises(ValueError):
#         OneTimeSigningKey(matrix1, matrix3)
#
# def test_OneTimeKeyTuple_class():
#     pass
#
#
# def test_SignatureChallenge_class():
#     pass
#
#
# def test_Signature_class():
#     pass
#
#
# def test_AggregationCoefficient_class():
#     pass
#
#
# def test_AggregateSignature_class():
#     pass
#
#
# # Test supporting internal functions.
#
# @pytest.mark.parametrize("secpar", [128, 256])
# def test_pre_hash_msg_to_int_digest(monkeypatch, example_params, example_msgs, secpar):
#     # Monkeypatch sha3_256_wrapper with dummy_sha3_256
#     monkeypatch.setattr('fusion.fusion.sha3_256_wrapper', dummy_sha3_256)
#
#     params = example_params[secpar]
#     msgs = example_msgs[secpar]
#
#     for msg in msgs:
#         hashed_int = _pre_hash_msg_to_int_digest(params, msg)
#         expected_int = int.from_bytes(b'\x00\x01' * 16, byteorder="little")
#         assert hashed_int == expected_int, f"Hashed int for message '{msg}' with secpar {secpar} was incorrect"
#
#
# def test_hash_vk_and_pre_hashed_msg_to_bytes_digest(mocker, example_params):
#     for secpar in [128, 256]:
#         # Create sample Params and OneTimeVerificationKey objects
#         params: Params = example_params[secpar]
#         otk: OneTimeKeyTuple = fusion_keygen(params)
#         otsk, otvk = otk.otsk, otk.otvk
#
#         pre_hashed_message = 1234567890
#
#         # Expected result, assuming the shake_256 function would return this value
#         expected_result = b"expected_shake_256_result"
#
#         # Create a mock object for the hash object returned by shake_256
#         mock_hash_obj = mocker.Mock()
#         mock_hash_obj.digest.return_value = expected_result
#
#         # Mock the shake_256 function to return the mock hash object
#         mock_shake_256 = mocker.patch(
#             "fusion.fusion.shake_256", return_value=mock_hash_obj
#         )
#
#         # Call the hash_vk_and_int_to_bytes function
#         # we are mocking a result, so n can be anything.
#         result = _hash_vk_and_pre_hashed_msg_to_bytes_digest(params=params, otvk=otvk, pre_hashed_msg=pre_hashed_message,
#                                                              num_bytes=0)
#
#         # Assert that the result matches the expected value
#         assert result == expected_result
#
#         # Check that the shake_256 function was called with the correct arguments
#         x = (
#             params.sign_hash_dst
#             + ","
#             + str(otvk)
#             + ","
#             + str(pre_hashed_message)
#         ).encode()
#         mock_shake_256.assert_called_once_with(x)
#
#
# @pytest.mark.parametrize("example_secpar, norm_bound, weight_bound, some_byte_value, expected_output", [
#     # Case 1: All zero bytes input for security parameter 128
#     (128, 1, 3, 0x00, [0, -1, -1,] + [0]*60 + [-1]),
#     # Case 2: All one bytes input for security parameter 128
#     (128, 2, 3, 0xFF, [0, -2, -2,] + [0]*60 + [-2]),
#     # Additional test cases as needed for different security parameters or bounds
# ])
# def test_bytes_digest_to_polys(example_params, example_secpar, norm_bound, weight_bound, some_byte_value,
#                                expected_output):
#     # Access the Params object corresponding to the security parameter example_secpar
#     params = example_params[example_secpar]
#     assert params.secpar == example_secpar
#
#     # Generate the input byte string based on 'one_byte_value'
#     total_bytes = bytes_for_bdd_poly(secpar=params.secpar, beta=norm_bound, degree=params.degree, omega=weight_bound)
#     input_bytes = bytes([some_byte_value]) * total_bytes
#
#     # Decode the input byte string
#     decoded_poly: List[Poly] = _bytes_digest_to_polys(params=params, beta=norm_bound, omega=weight_bound, b=input_bytes, num_polys=1)
#
#     # Check if the decoded coefficients match the expected output
#     assert decoded_poly[0].coefficients == expected_output, f"Expected {expected_output}, but got {decoded_poly}"
#
#
# def test_make_sig_challenge():
#     pass
#
#
# def test_verify_signature():
#     pass
#
#
# def test_hash_vks_and_pre_hashed_msgs_and_challs_to_bytes_digest():
#     pass
#
#
# def test_make_aggregation_coefficients():
#     pass
#
#
# # Fixtures
#
# @pytest.fixture
# def example_params() -> Dict[int, Params]:
#     # Create and return a params object
#     return {128: Params(secpar=128), 256: Params(secpar=256)}
#
#
# @pytest.fixture
# def example_keys(example_params) -> Dict[int, List[OneTimeKeyTuple]]:
#     return {128: [fusion_keygen(example_params[128]) for _ in range(3)], 256: [fusion_keygen(example_params[256]) for _ in range(3)]}
#
#
# @pytest.fixture
# def example_msgs(example_params, example_keys) -> Dict[int, List[bytes]]:
#     return {128: [f"example message {i} for 128 bit signing" for i in range(3)], 256: [f"example message {i} for 256 bit signing" for i in range(3)]}
#
#
# @pytest.fixture
# def example_challenges(example_params, example_keys, example_msgs) -> Dict[int, List[SignatureChallenge]]:
#     return {128: [_make_sig_challenge(params=example_params[128], otvk=otk.otvk, msg=msg) for otk, msg in zip(example_keys[128], example_msgs[128])],
#             256: [_make_sig_challenge(params=example_params[256], otvk=otk.otvk, msg=msg) for otk, msg in zip(example_keys[256], example_msgs[256])]}
#
#
# # Test fusion algorithms
#
# def test_fusion_setup():
#     def assert_params_equality_except_public_challenges(params_instance_1, params_instance_2):
#         for attr in vars(params_instance_1):
#             if attr == 'public_challenge':
#                 continue  # Skip this attribute
#             attr_value_1 = getattr(params_instance_1, attr)
#             attr_value_2 = getattr(params_instance_2, attr)
#             assert attr_value_1 == attr_value_2
#
#     for next_secpar in [128, 256]:
#         expected_params: Params = Params(secpar=next_secpar)
#         observed_params: Params = fusion_setup(secpar=next_secpar)
#         assert_params_equality_except_public_challenges(expected_params, observed_params)
#
#
# def test_fusion_keygen(example_params):
#     for secpar in [128, 256]:
#         # Setup
#         params: Params = example_params[secpar]
#
#         # Keygen
#         otk: OneTimeKeyTuple = fusion_keygen(params=params)
#         assert isinstance(otk, OneTimeKeyTuple)
#         assert hasattr(otk, "otsk")
#         assert hasattr(otk, "otvk")
#         assert isinstance(otk.otsk, OneTimeSigningKey)
#         assert isinstance(otk.otvk, OneTimeVerificationKey)
#
#         # Unpack the key
#         otsk: OneTimeSigningKey  # private
#         otvk: OneTimeVerificationKey  # public
#         otsk, otvk = otk.otsk, otk.otvk  # unpack
#
#         assert hasattr(otsk, "left_sk_hat")
#         assert hasattr(otsk, "right_sk_hat")
#         assert isinstance(otsk.left_sk_hat, GeneralMatrix)
#         assert isinstance(otsk.right_sk_hat, GeneralMatrix)
#
#         inv_left_sk_hat: GeneralMatrix = GeneralMatrix(
#             matrix=[[transform(f) for f in row] for row in otsk.left_sk_hat.matrix]
#         )
#         inv_right_sk_hat: GeneralMatrix = GeneralMatrix(
#             matrix=[[transform(f) for f in row] for row in otsk.right_sk_hat.matrix]
#         )
#         assert inv_left_sk_hat.norm(p="infty") <= params.beta_sk
#         assert inv_right_sk_hat.norm(p="infty") <= params.beta_sk
#         assert inv_left_sk_hat.weight() <= params.omega_sk
#         assert inv_right_sk_hat.weight() <= params.omega_sk
#
#         assert hasattr(otvk, "left_vk_hat")
#         assert hasattr(otvk, "right_vk_hat")
#         assert isinstance(otvk.left_vk_hat, GeneralMatrix)
#         assert isinstance(otvk.right_vk_hat, GeneralMatrix)
#
#         assert params.public_challenge * otsk.left_sk_hat == otvk.left_vk_hat
#         assert params.public_challenge * otsk.right_sk_hat == otvk.right_vk_hat
#
#
# def test_fusion_sign():
#     pass
#
#
# def test_fusion_aggregate():
#     pass
#
#
# def test_fusion_verify():
#     pass