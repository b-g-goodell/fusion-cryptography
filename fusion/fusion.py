# from hashlib import shake_256, sha3_256
# from math import ceil, log2
# from typing import List, Dict, Tuple, Union
# from algebra.matrices import GeneralMatrix
# from algebra.polynomials import PolynomialCoefficientRepresentation, PolynomialNTTRepresentation, transform
# from algebra.sampling import sample_coefficient_matrix, sample_ntt_matrix
#
# # Some global constants
#
# PRIME: int = 2147465729
# DEGREES: Dict[int, int] = {128: 2**6, 256: 2**8}
# RANKS: Dict[int, int] = {128: 195, 256: 83}
# CAPACITIES: Dict[int, int] = {128: 1796, 256: 2818}
# ROOTS: Dict[int, int] = {128: 23584283, 256: 3337519}
# SIGN_PRE_HASH_DSTS: Dict[int, str] = {128: '0000', 256: '0001'}
# SIGN_HASH_DSTS: Dict[int, str] = {128: '0010', 256: '0011'}
# AGG_XOF_DSTS: Dict[int, str] = {128: '0100', 256: '0101'}
# CH_WTS: Dict[int, int] = {128: 27, 256: 60}
# AG_WTS: Dict[int, int] = {128: 35, 256: 60}
# CH_BDS: Dict[int, int] = {128: 3, 256: 1}
# AG_BDS: Dict[int, int] = {128: 2, 256: 1}
# PREFIX_PARAMETERS: Dict[int, Dict[str, Union[int, str]]] = {secpar: {
#     "capacity": CAPACITIES[secpar], "modulus": PRIME, "degree": DEGREES[secpar],
#     "root_order": 2 * DEGREES[secpar], "root": ROOTS[secpar], "inv_root": pow(ROOTS[secpar], PRIME - 2, PRIME),
#     "num_rows_pub_challenge": 1, "num_rows_sk": RANKS[secpar], "num_rows_vk": 1,
#     "num_cols_sk": 1, "num_cols_vk": 1, "sign_pre_hash_dst": SIGN_PRE_HASH_DSTS[secpar],
#     "sign_hash_dst": SIGN_HASH_DSTS[secpar], "agg_xof_dst": AGG_XOF_DSTS[secpar], "beta_sk": 52, "beta_ch": 1,
#     "beta_ag": 1, "omega_sk": DEGREES[secpar], "omega_ch": CH_WTS[secpar], "omega_ag": AG_WTS[secpar]
# } for secpar in [128, 256]}
#
# # Constants computed from other constants
#
# for secpar in [128, 256]:
#     PREFIX_PARAMETERS[secpar]["omega_vf_intermediate"] = max(0, min(PREFIX_PARAMETERS[secpar]["degree"], PREFIX_PARAMETERS[secpar]["omega_sk"]*(1+PREFIX_PARAMETERS[secpar]["omega_ch"])))
#     PREFIX_PARAMETERS[secpar]["beta_vf_intermediate"] = PREFIX_PARAMETERS[secpar]["beta_sk"] * (1 + max(0, min(PREFIX_PARAMETERS[secpar]["degree"], PREFIX_PARAMETERS[secpar]["omega_ch"]))*PREFIX_PARAMETERS[secpar]["beta_ch"])
#
#     PREFIX_PARAMETERS[secpar]["omega_vf"] = max(0, min(PREFIX_PARAMETERS[secpar]["degree"], PREFIX_PARAMETERS[secpar]["capacity"]*PREFIX_PARAMETERS[secpar]["omega_vf_intermediate"]))
#     PREFIX_PARAMETERS[secpar]["beta_vf"] = (PREFIX_PARAMETERS[secpar]["capacity"] * max(0, min(PREFIX_PARAMETERS[secpar]["degree"], PREFIX_PARAMETERS[secpar]["omega_ag"])) * PREFIX_PARAMETERS[secpar]["beta_ag"] * PREFIX_PARAMETERS[secpar]["beta_vf_intermediate"])
#
# # Simple wrappers for sha3 and shake
#
# def sha3_256_wrapper(message: str | bytes) -> bytes:
#     """ Simple wrapper for sha3_256 """
#     if isinstance(message, str):
#         return sha3_256(message.encode('utf-8')).digest()
#     return sha3_256(message).digest()
#
#
# def shake_256_wrapper(x: str | bytes, n: int) -> bytes:
#     """ Simple wrapper for shake_256 """
#     if isinstance(x, str):
#         return shake_256(x.encode('utf-8')).digest(n)
#     return shake_256(x).digest(n)
#
#
# # simple functions for counting bits and bytes
#
# def bits_for_bdd_coef(secpar: int, beta: int) -> int:
#     """
#     Number of bits to sample an integer from list(range(2*beta+1)) with
#     bias (statistical distance from uniformity) which is O(2**-secpar)
#     """
#     return ceil(log2(2*beta+1)) + secpar
#
#
# def bytes_for_bdd_coef(secpar: int, beta: int) -> int:
#     """
#     Number of bytes to sample an integer from list(range(2*beta+1)) with
#     bias (statistical distance from uniformity) which is O(2**-secpar)
#     """
#     return ceil(bits_for_bdd_coef(secpar=secpar, beta=beta)/8)
#
#
# def bits_for_index(secpar: int, degree: int) -> int:
#     """
#     Number of bits to sample a monomial degree from list(range(degree))
#     with bias (statistical distance from uniformity) which is
#     O(2**-secpar)
#     """
#     if degree < 2:
#         raise ValueError("Must have degree >= 2")
#     return ceil(log2(degree)) + secpar
#
#
# def bytes_for_index(secpar: int, degree: int) -> int:
#     """
#     Number of bytes to sample a monomial X**j with exponent j from
#     list(range(degree)) with bias (statistical distance from uni-
#     formity) which is O(2**-secpar)
#     """
#     return ceil(bits_for_index(secpar=secpar, degree=degree)/8)
#
#
# def bits_for_fy_shuffle(secpar: int, degree: int) -> int:
#     """
#     Number of bits to Fisher-Yates shuffle list(range(degree)) with
#     bias (statistical distance from uniformity) which is O(2**-secpar)
#     """
#     return degree * bits_for_index(secpar=secpar, degree=degree)
#
#
# def bytes_for_fy_shuffle(secpar: int, degree: int) -> int:
#     """
#     Number of bytes to Fisher-Yates shuffle list(range(degree)) with
#     bias (statistical distance from uniformity) which is O(2**-secpar)
#     """
#     return ceil(bits_for_fy_shuffle(secpar=secpar, degree=degree)/8)
#
#
# def bits_for_bdd_poly(secpar: int, beta: int, degree: int, omega: int) -> int:
#     """
#     Number of bits to sample a polynomial which is a sum of omega
#     monomials whose coefficients are in list(range(-beta,beta+1))
#     with bias (statistical distance from uniformity) which is
#     O(2**-secpar)
#     """
#     return omega * bits_for_bdd_coef(secpar=secpar, beta=beta) + bits_for_fy_shuffle(secpar=secpar, degree=degree)
#
#
# def bytes_for_bdd_poly(secpar: int, beta: int, degree: int, omega: int) -> int:
#     """
#     Number of bytes to sample a polynomial which is a sum of omega
#     monomials whose coefficients are in list(range(-beta,beta+1))
#     with bias (statistical distance from uniformity) which is
#     O(2**-secpar)
#     """
#     return ceil(bits_for_bdd_poly(secpar=secpar, beta=beta, degree=degree, omega=omega) / 8)
#
# # Typing
#
# class Params(object):
#     secpar: int
#     capacity: int
#     modulus: int
#     degree: int
#     root_order: int
#     root: int
#     inv_root: int
#     num_rows_pub_challenge: int
#     num_rows_sk: int
#     num_rows_vk: int
#     num_cols_pub_challenge: int
#     num_cols_sk: int
#     num_cols_vk: int
#     beta_sk: int
#     beta_ch: int
#     beta_ag: int
#     beta_vf_intermediate: int
#     beta_vf: int
#     omega_sk: int
#     omega_ch: int
#     omega_ag: int
#     omega_vf_intermediate: int
#     omega_vf: int
#     public_challenge: GeneralMatrix
#     sign_pre_hash_dst: str
#     sign_hash_dst: str
#     agg_xof_dst: str
#     bytes_for_one_coef_bdd_by_beta_ch: int
#     bytes_for_one_coef_bdd_by_beta_ag: int
#     bytes_for_fy_shuffle: int
#     bytes_per_sig_chall_poly: int
#     bytes_per_agg_coef_poly: int
#
#     def __init__(self, secpar: int):
#         if secpar not in PREFIX_PARAMETERS:
#             raise ValueError("Invalid security parameter.")
#         self.secpar = secpar
#         self.capacity = PREFIX_PARAMETERS[secpar]["capacity"]
#         self.modulus = PREFIX_PARAMETERS[secpar]["modulus"]
#         self.degree = PREFIX_PARAMETERS[secpar]["degree"]
#         self.root_order = PREFIX_PARAMETERS[secpar]["root_order"]
#         self.root = PREFIX_PARAMETERS[secpar]["root"]
#         self.inv_root = PREFIX_PARAMETERS[secpar]["inv_root"]
#         self.num_rows_pub_challenge = PREFIX_PARAMETERS[secpar]["num_rows_pub_challenge"]
#         self.num_rows_sk = PREFIX_PARAMETERS[secpar]["num_rows_sk"]
#         self.num_rows_vk = PREFIX_PARAMETERS[secpar]["num_rows_vk"]
#         self.num_cols_sk = PREFIX_PARAMETERS[secpar]["num_cols_sk"]
#         self.num_cols_vk = PREFIX_PARAMETERS[secpar]["num_cols_vk"]
#         self.sign_pre_hash_dst = PREFIX_PARAMETERS[secpar]["sign_pre_hash_dst"]
#         self.sign_hash_dst = PREFIX_PARAMETERS[secpar]["sign_hash_dst"]
#         self.agg_xof_dst = PREFIX_PARAMETERS[secpar]["agg_xof_dst"]
#         self.beta_sk = PREFIX_PARAMETERS[secpar]["beta_sk"]
#         self.beta_ch = PREFIX_PARAMETERS[secpar]["beta_ch"]
#         self.beta_ag = PREFIX_PARAMETERS[secpar]["beta_ag"]
#         self.beta_vf_intermediate = PREFIX_PARAMETERS[secpar]["beta_vf_intermediate"]
#         self.beta_vf = PREFIX_PARAMETERS[secpar]["beta_vf"]
#         self.omega_sk = PREFIX_PARAMETERS[secpar]["omega_sk"]
#         self.omega_ch = PREFIX_PARAMETERS[secpar]["omega_ch"]
#         self.omega_ag = PREFIX_PARAMETERS[secpar]["omega_ag"]
#         self.omega_vf_intermediate = PREFIX_PARAMETERS[secpar]["omega_vf_intermediate"]
#         self.omega_vf = PREFIX_PARAMETERS[secpar]["omega_vf"]
#         self.public_challenge = sample_ntt_matrix(modulus=self.modulus, degree=self.degree, root_order=self.root_order, root=self.root, inv_root=self.inv_root, num_rows=self.num_rows_pub_challenge, num_cols=self.num_rows_sk)
#
#     def __str__(self) -> str:
#         return f"Params(secpar={self.secpar}, capacity={self.capacity}, modulus={self.modulus}, degree={self.degree}, root_order={self.root_order}, root={self.root}, inv_root={self.inv_root}, num_rows_pub_challenge={self.num_rows_pub_challenge}, num_rows_sk={self.num_rows_sk}, num_rows_vk={self.num_rows_vk}, num_cols_sk={self.num_cols_sk}, num_cols_vk={self.num_cols_vk}, beta_sk={self.beta_sk}, beta_ch={self.beta_ch}, beta_ag={self.beta_ag}, beta_vf={self.beta_vf}, omega_sk={self.omega_sk}, omega_ch={self.omega_ch}, omega_ag={self.omega_ag}, omega_vf={self.omega_vf}, public_challenge={str(self.public_challenge)}, sign_pre_hash_dst={self.sign_pre_hash_dst}, sign_hash_dst={self.sign_hash_dst}, agg_xof_dst={self.agg_xof_dst}, bytes_for_one_coef_bdd_by_beta_ch={self.bytes_for_one_coef_bdd_by_beta_ch}, bytes_for_one_coef_bdd_by_beta_ag={self.bytes_for_one_coef_bdd_by_beta_ag}, bytes_for_fy_shuffle={self.bytes_for_fy_shuffle})"
#
#     def __repr__(self) -> str:
#         return self.__str__()
#
#     def __eq__(self, other):
#         return self.__dict__ == other.__dict__
#
#
# class OneTimeSigningKey(object):
#     # Todo: encoding/decoding separate from sampling
#     left_sk_hat: GeneralMatrix
#     right_sk_hat: GeneralMatrix
#
#     def __init__(
#         self,
#         left_sk_hat: GeneralMatrix,
#         right_sk_hat: GeneralMatrix,
#     ):
#         if not isinstance(left_sk_hat, GeneralMatrix):
#             raise TypeError("left_sk_hat must be GeneralMatrix")
#         elif not isinstance(right_sk_hat, GeneralMatrix):
#             raise TypeError("right_sk_hat must be GeneralMatrix")
#         elif left_sk_hat.elem_class != right_sk_hat.elem_class:
#             raise ValueError("Both GeneralMatrix inputs must have same elementary class.")
#         elif len(left_sk_hat.matrix) != len(right_sk_hat.matrix) or len(left_sk_hat.matrix[0]) != len(right_sk_hat.matrix[0]):
#             raise ValueError("Dimension mismatch.")
#         self.left_sk_hat = left_sk_hat
#         self.right_sk_hat = right_sk_hat
#
#     def __str__(self):
#         # Do not print your key.
#         return f"OneTimeSigningKey(left_sk_hat={str(self.left_sk_hat)}, right_sk_hat={str(self.right_sk_hat)})"
#
#     def __repr__(self):
#         # Do not print your key.
#         return self.__str__()
#
#
# class OneTimeVerificationKey(object):
#     # Todo: encoding/decoding separate from sampling
#     left_vk_hat: GeneralMatrix
#     right_vk_hat: GeneralMatrix
#
#     def __init__(self, left_vk_hat: GeneralMatrix, right_vk_hat: GeneralMatrix):
#         if not isinstance(left_vk_hat, GeneralMatrix):
#             raise TypeError("left_sk_hat must be GeneralMatrix")
#         elif not isinstance(right_vk_hat, GeneralMatrix):
#             raise TypeError("right_sk_hat must be GeneralMatrix")
#         elif left_vk_hat.elem_class != right_vk_hat.elem_class:
#             raise ValueError("Both GeneralMatrix inputs must have same elementary class.")
#         elif len(left_vk_hat.matrix) != len(right_vk_hat.matrix) or len(left_vk_hat.matrix[0]) != len(right_vk_hat.matrix[0]):
#             raise ValueError("Dimension mismatch.")
#         self.left_vk_hat = left_vk_hat
#         self.right_vk_hat = right_vk_hat
#
#     def __str__(self):
#         return f"OneTimeVerificationKey(left_vk_hat={self.left_vk_hat}, right_vk_hat={self.right_vk_hat})"
#
#     def __repr__(self):
#         return self.__str__()
#
#
# class OneTimeKeyTuple(object):
#     otsk: OneTimeSigningKey
#     otvk: OneTimeVerificationKey
#
#     def __init__(self, otsk: OneTimeSigningKey, otvk: OneTimeVerificationKey):
#         self.otsk = otsk
#         self.otvk = otvk
#
#     def __str__(self):
#         return f"OneTimeKeyTuple(otsk={self.otsk},otvk={self.otvk})"
#
#     def __repr__(self):
#         return self.__str__()
#
#
# class SignatureChallenge(object):
#     # Todo: encoding/decoding separate from sampling
#     c_hat: PolynomialNTTRepresentation
#
#     def __init__(self, c_hat: PolynomialNTTRepresentation):
#         if not isinstance(c_hat, PolynomialNTTRepresentation):
#             raise TypeError(f"c_hat must be a PolynomialNTTRepresentation object but had type(c_hat)={type(c_hat)}")
#         self.c_hat = c_hat
#
#     def __str__(self):
#         return f"SignatureChallenge(c_hat={str(self.c_hat)})"
#
#     def __repr__(self):
#         return self.__str__()
#
#     def __eq__(self, other):
#         return self.c_hat == other.c_hat
#
#
# class Signature(object):
#     signature_hat: GeneralMatrix
#
#     def __init__(self, signature_hat: GeneralMatrix):
#         self.signature_hat = signature_hat
#
#     def __str__(self):
#         return f"Signature(signature_hat={str(self.signature_hat)})"
#
#     def __repr__(self):
#         return self.__str__()
#
#
# class AggregationCoefficient(object):
#     alpha_hat: PolynomialNTTRepresentation
#
#     def __init__(self, alpha_hat: PolynomialNTTRepresentation):
#         self.alpha_hat = alpha_hat
#
#     def __str__(self):
#         return f"AggregationCoefficient(alpha_hat={self.alpha_hat})"
#
#     def __repr__(self):
#         return self.__str__()
#
#
# class AggregateSignature(object):
#     aggregate_signature_hat: GeneralMatrix
#
#     def __init__(self, aggregate_signature_hat: GeneralMatrix):
#         self.aggregate_signature_hat = aggregate_signature_hat
#
#     def __str__(self):
#         return f"AggregateSignature(aggregate_signature_hat={str(self.aggregate_signature_hat)})"
#
#     def __repr__(self):
#         return self.__str__()
#
#
# # Supporting internal functions
#
# def _pre_hash_msg_to_int_digest(params: Params, msg: str) -> int:
#     # Apply sha3-256 for hash-then-sign
#     dst_and_msg: bytes = str.encode(
#         params.sign_pre_hash_dst
#         + ","
#         + msg
#     )
#     pre_hashed_message: bytes = sha3_256_wrapper(dst_and_msg)
#     pre_hashed_message_as_int: int = int.from_bytes(pre_hashed_message, byteorder="little")
#     return pre_hashed_message_as_int
#
#
# def _hash_vk_and_pre_hashed_msg_to_bytes_digest(params: Params, otvk: OneTimeVerificationKey, pre_hashed_msg: int, num_bytes: int) -> bytes:
#     dst_and_otk_and_prehashed_msg: bytes = str.encode(
#         params.sign_hash_dst
#         + ","
#         + str(otvk)
#         + ","
#         + str(pre_hashed_msg)
#     )
#     return shake_256_wrapper(dst_and_otk_and_prehashed_msg, num_bytes)
#
#
# def _bytes_digest_to_polys(params: Params, beta: int, omega: int, b: bytes, num_polys: int) -> List[PolynomialCoefficientRepresentation]:
#     if omega < 1:
#         raise ValueError("Cannot decode bytes to no coefficients.")
#     elif beta < 1:
#         raise ValueError("Cannot decode bytes to non-positive coefficients.")
#     elif num_polys < 1:
#         raise ValueError("Cannot decode to no polynomials.")
#     coef_bytes: int = bytes_for_bdd_coef(secpar=params.secpar, beta=beta)
#     idx_bytes: int = bytes_for_index(secpar=params.secpar, degree=params.degree)
#     bytes_per_bdd_poly: int = bytes_for_bdd_poly(secpar=params.secpar, beta=beta, degree=params.degree, omega=omega)
#     if len(b) < num_polys * bytes_per_bdd_poly:
#         raise ValueError(
#             f"Too few bytes to decode polynomial. Expected {bytes_per_bdd_poly} but got {len(b)}"
#         )
#     result: List[PolynomialCoefficientRepresentation] = []
#     next_b: bytes
#     remaining_b: bytes = b
#     while remaining_b or next_b:
#         next_b, remaining_b = (
#             remaining_b[:bytes_per_bdd_poly],
#             remaining_b[bytes_per_bdd_poly:]
#         )
#         # Generate the coefficients
#         coefficients: List[int] = []
#         next_bytes: bytes
#         remaining_bytes: bytes = remaining_b
#         for i in range(omega):
#             next_bytes, remaining_bytes = (
#                 remaining_bytes[:coef_bytes],
#                 remaining_bytes[coef_bytes:],
#             )
#             coefficients += [
#                 (int.from_bytes(next_bytes, byteorder="little") % (2*beta + 1)) - beta
#             ]
#         while len(coefficients) < params.degree:
#             coefficients += [0]
#
#         # Fisher-Yates shuffle if necessary
#         # This part can be optimized a lot
#         if omega < params.degree:
#             for i in range(params.degree - 1, omega, -1):
#                 # convert next_few_bytes to an integer modulo i+1
#                 next_bytes, remaining_bytes = (
#                     remaining_bytes[:idx_bytes],
#                     remaining_bytes[idx_bytes:],
#                 )
#                 j: int = int.from_bytes(next_bytes, byteorder="little") % (i + 1)
#                 coefficients[i], coefficients[j] = coefficients[j], coefficients[i]
#         result += [PolynomialCoefficientRepresentation(modulus=params.modulus, degree=params.degree, root=params.root,
#                                                inv_root=params.inv_root, root_order=params.root_order,
#                                                coefficients=coefficients)]
#     return result
#
#
# def _make_sig_challenge(params: Params, otvk: OneTimeVerificationKey, msg: str) -> SignatureChallenge:
#     num_bytes_needed: int = bytes_for_bdd_poly((params.secpar, params.beta_ch, params.degree, params.omega_ch))
#     prehashed_msg: int = _pre_hash_msg_to_int_digest(params=params, msg=msg)
#     hashed_key_and_prehashed_msg: bytes = _hash_vk_and_pre_hashed_msg_to_bytes_digest(params=params, otvk=otvk,
#                                                                                       pre_hashed_msg=prehashed_msg,
#                                                                                       num_bytes=num_bytes_needed)
#     c: PolynomialCoefficientRepresentation = _bytes_digest_to_polys(params=params, beta=params.beta_ch, omega=params.omega_ch, b=hashed_key_and_prehashed_msg, num_polys=1)
#     return SignatureChallenge(c_hat=transform(c))
#
#
# def _verify_signature(params: Params, otvk: OneTimeVerificationKey, msg: str, sig: Signature) -> Tuple[bool, str]:
#     if sig.signature_hat.norm() > params.beta_vf_intermediate:
#         return False, "Norm of purported signature too large"
#     elif sig.signature_hat.weight() > params.omega_vf_intermediate:
#         return False, "Hamming weight of purported signature too large"
#     c_hat: SignatureChallenge = _make_sig_challenge(params=params, otvk=otvk, msg=msg)
#     target_image: GeneralMatrix = otvk.left_vk_hat + c_hat * otvk.right_vk_hat
#     if params.public_challenge * sig != target_image:
#         return False, "Signature image does not match target image."
#     return True, "Signature valid."
#
#
# def _hash_vks_and_pre_hashed_msgs_and_challs_to_bytes_digest(params: Params, otvks: List[OneTimeVerificationKey],
#                                                              pre_hashed_msgs: List[int], challs: List[SignatureChallenge],
#                                                              num_bytes: int) -> bytes:
#     dst_and_vks_and_pre_hashed_msgs_and_challs: bytes = str.encode(
#         params.agg_xof_dst
#         + ","
#         + str(list(zip(otvks, pre_hashed_msgs, challs)))
#     )
#     return shake_256_wrapper(dst_and_vks_and_pre_hashed_msgs_and_challs, num_bytes)
#
#
# def _make_aggregation_coefficients(params: Params, srt_otvks: List[OneTimeVerificationKey],
#                                    srt_pre_hashed_msgs: List[int], srt_challs: List[SignatureChallenge]) -> List[AggregationCoefficient]:
#     num_keys: int = len(srt_otvks)
#     num_bytes_needed_per_poly: int = bytes_for_bdd_poly((params.secpar, params.beta_ag, params.degree, params.omega_ag))
#     num_bytes_needed: int = num_keys * num_bytes_needed_per_poly
#     bytes_digest: bytes = _hash_vks_and_pre_hashed_msgs_and_challs_to_bytes_digest(params=params, otvks=srt_otvks, pre_hashed_msgs=srt_pre_hashed_msgs, challs=srt_challs, num_bytes=num_bytes_needed)
#     alphas: List[PolynomialCoefficientRepresentation] = _bytes_digest_to_polys(params=params, beta=params.beta_ag, degree=params.degree, omega=params.omega_ag, b=bytes_digest, num_polys=num_keys)
#     alpha_hats: List[PolynomialNTTRepresentation] = [transform(a) for a in alphas]
#     return [AggregationCoefficient(alpha_hat=a_hat) for a_hat in alpha_hats]
#
#
# # Fusion algorithms
#
# def fusion_setup(secpar: int) -> Params:
#     return Params(secpar=secpar)
#
#
# def fusion_keygen(params: Params) -> OneTimeKeyTuple:
#     left_key_coefs: GeneralMatrix = sample_coefficient_matrix(modulus=params.modulus, degree=params.degree,
#                                                               root_order=params.root_order, root=params.root,
#                                                               inv_root=params.inv_root, num_rows=params.num_rows_sk,
#                                                               num_cols=params.num_cols_sk, norm_bound=params.beta_sk,
#                                                               weight_bound=params.omega_sk)
#     right_key_coefs: GeneralMatrix = sample_coefficient_matrix(modulus=params.modulus, degree=params.degree,
#                                                                root_order=params.root_order, root=params.root,
#                                                                inv_root=params.inv_root, num_rows=params.num_rows_sk,
#                                                                num_cols=params.num_cols_sk, norm_bound=params.beta_sk,
#                                                                weight_bound=params.omega_sk)
#     left_sk_hat: GeneralMatrix = GeneralMatrix(
#         matrix=[[transform(y) for y in z] for z in left_key_coefs.matrix]
#     )
#     right_sk_hat: GeneralMatrix = GeneralMatrix(
#         matrix=[[transform(y) for y in z] for z in right_key_coefs.matrix]
#     )
#     otsk: OneTimeSigningKey = OneTimeSigningKey(left_sk_hat=left_sk_hat, right_sk_hat=right_sk_hat)
#     left_vk_hat: GeneralMatrix = params.public_challenge * left_sk_hat
#     right_vk_hat: GeneralMatrix = params.public_challenge * right_sk_hat
#     otvk: OneTimeVerificationKey = OneTimeVerificationKey(left_vk_hat=left_vk_hat, right_vk_hat=right_vk_hat)
#     return OneTimeKeyTuple(otsk=otsk, otvk=otvk)
#
#
# def fusion_sign(params: Params, otk: OneTimeKeyTuple, msg: str) -> Signature:
#     c_hat: SignatureChallenge = _make_sig_challenge(params=params, otvk=otk.otvk, msg=msg)
#     return Signature(signature_hat=otk.otsk.left_sk_hat * c_hat.c_hat + otk.otsk.right_sk_hat)
#
#
# # def fusion_aggregate(params: Params, otvks: List[OneTimeVerificationKey], msgs: List[str], sigs: List[Signature]) -> AggregateSignature:
# #     pre_hashed_messages: List[int] = [_pre_hash_msg_to_int_digest(params=params, msg=msg) for msg in msgs]
# #     challs: List[SignatureChallenge] = [_make_sig_challenge(params=params, otvk=otvk, msg=msg) for otvk, msg in zip(otvks, msgs)]
# #     sorted_input: List[Tuple[OneTimeVerificationKey, str, Signature, int, SignatureChallenge]] = sorted(
# #         list(zip(otvks, msgs, sigs, pre_hashed_messages, challs)), key=lambda x: str(x[0])
# #     )
# #     sorted_otvks: List[Signature] = [x[0] for x in sorted_input]
# #     sorted_signatures: List[Signature] = [x[2] for x in sorted_input]
# #     sorted_prehashed_messages: List[int] = [x[3] for x in sorted_input]
# #     sorted_challs: List[SignatureChallenge] = [x[4] for x in sorted_input]
# #
# #     aggregation_coefficients: List[AggregationCoefficient] = _make_aggregation_coefficients(params=params,
# #                                                                                             srt_otvks=sorted_otvks,
# #                                                                                             srt_pre_hashed_msgs=sorted_prehashed_messages,
# #                                                                                             srt_challs=sorted_challs)
# #     aggregate_signature_hat_values: GeneralMatrix = (sorted_signatures[0].signature_hat * aggregation_coefficients[0].alpha_hat)
# #     for next_alpha, next_sig in zip(aggregation_coefficients[1:], sorted_signatures[1:]):
# #         aggregate_signature_hat_values += next_sig.signature_hat * next_alpha.alpha_hat
# #     return Signature(signature_hat=aggregate_signature_hat_values)
#
# #
# # def fusion_verify(
# #     params: Params,
# #     keys: List[OneTimeVerificationKey],
# #     messages: List[str],
# #     aggregate_signature: Signature,
# # ) -> Tuple[bool, str]:
# #     if len(keys) > params.capacity:
# #         return False, f"Too many keys."
# #     elif len(keys) != len(messages):
# #         return False, f"Number of keys and messages must be equal."
# #     coef_rep_agg_sig: GeneralMatrix = GeneralMatrix(
# #         matrix=[[transform(z) for z in y] for y in aggregate_signature.signature_hat]
# #     )
# #     sorted_input = sorted(zip(keys, messages), key=lambda x: str(x[0]))
# #     sorted_vks: List[OneTimeVerificationKey] = [x[0] for x in sorted_input]
# #     sorted_challs: List[SignatureChallenge] = [
# #         _hash_ch(params=params, key=next_vk, message=next_m)
# #         for next_vk, next_m in sorted_input
# #     ]
# #     aggregation_coefficients: List[AggregationCoefficient] = _make_agg_coefs(params=params,
# #                                                                              keys=[x[0] for x in sorted_input],
# #                                                                              messages=[x[1] for x in sorted_input])
# #     tmp: GeneralMatrix = sorted_vks[0].left_vk_hat * sorted_challs[0].c_hat
# #     tmp += sorted_vks[0].right_vk_hat
# #     target: GeneralMatrix = (
# #         sorted_vks[0].left_vk_hat * sorted_challs[0].c_hat + sorted_vks[0].right_vk_hat
# #     ) * aggregation_coefficients[0].alpha_hat
# #     for next_alpha, next_vk, next_chall in zip(
# #         aggregation_coefficients[1:], sorted_vks[1:], sorted_challs[1:]
# #     ):
# #         target += (
# #             next_vk.left_vk_hat * next_chall.c_hat + next_vk.right_vk_hat
# #         ) * next_alpha.alpha_hat
# #     observed: GeneralMatrix = (
# #         params.public_challenge * aggregate_signature.signature_hat
# #     )
# #     for a, b in zip(target.matrix, observed.matrix):
# #         for c, d in zip(a, b):
# #             if c != d:
# #                 return False, f"Target doesn't match image of aggregate signature."
# #     if any(
# #         z.norm(p="infty") > params.beta_vf for y in coef_rep_agg_sig.matrix for z in y
# #     ):
# #         return False, f"Norm of aggregate signature too large."
# #     elif any(z.weight() > params.omega_vf for y in coef_rep_agg_sig.matrix for z in y):
# #         return False, f"Weight of aggregate signature too large."
# #     return True, ""
