from hashlib import shake_256, sha3_256
from math import ceil, log2
from algebra.polynomials import _Polynomial as Poly
from algebra.matrices import _PolynomialMatrix as Mat
from crypto.fusion.params import IACR_SUGGESTED_PARAMS


# Constants computed from other constants
for secpar in [128, 256]:
    IACR_SUGGESTED_PARAMS[secpar]["omega_vf_intermediate"] = max(0, min(IACR_SUGGESTED_PARAMS[secpar]["degree"], IACR_SUGGESTED_PARAMS[secpar]["omega_sk"] * (1 + IACR_SUGGESTED_PARAMS[secpar]["omega_ch"])))
    IACR_SUGGESTED_PARAMS[secpar]["beta_vf_intermediate"] = IACR_SUGGESTED_PARAMS[secpar]["beta_sk"] * (1 + max(0, min(IACR_SUGGESTED_PARAMS[secpar]["degree"], IACR_SUGGESTED_PARAMS[secpar]["omega_ch"])) * IACR_SUGGESTED_PARAMS[secpar]["beta_ch"])

    IACR_SUGGESTED_PARAMS[secpar]["omega_vf"] = max(0, min(IACR_SUGGESTED_PARAMS[secpar]["degree"], IACR_SUGGESTED_PARAMS[secpar]["capacity"] * IACR_SUGGESTED_PARAMS[secpar]["omega_vf_intermediate"]))
    IACR_SUGGESTED_PARAMS[secpar]["beta_vf"] = (IACR_SUGGESTED_PARAMS[secpar]["capacity"] * max(0, min(IACR_SUGGESTED_PARAMS[secpar]["degree"], IACR_SUGGESTED_PARAMS[secpar]["omega_ag"])) * IACR_SUGGESTED_PARAMS[secpar]["beta_ag"] * IACR_SUGGESTED_PARAMS[secpar]["beta_vf_intermediate"])


# Simple wrappers for sha3 and shake
def sha3_256_wrapper(message: str | bytes) -> bytes:
    """ Simple wrapper for sha3_256 """
    if isinstance(message, str):
        return sha3_256(message.encode('utf-8')).digest()
    return sha3_256(message).digest()


def shake_256_wrapper(message: str | bytes, num_bytes: int) -> bytes:
    """ Simple wrapper for shake_256 """
    if isinstance(message, str):
        return shake_256(message.encode('utf-8')).digest(num_bytes)
    return shake_256(message).digest(num_bytes)


# counting bits and bytes
def bits_for_bdd_coef(secpar: int, beta: int) -> int:
    """
    Number of bits to sample an integer from list(range(2*beta+1)) with
    bias (statistical distance from uniformity) which is O(2**-secpar)
    """
    return ceil(log2(2*beta+1)) + secpar


def bytes_for_bdd_coef(secpar: int, beta: int) -> int:
    """
    Number of bytes to sample an integer from list(range(2*beta+1)) with
    bias (statistical distance from uniformity) which is O(2**-secpar)
    """
    return ceil(bits_for_bdd_coef(secpar=secpar, beta=beta)/8)


def bits_for_index(secpar: int, degree: int) -> int:
    """
    Number of bits to sample a monomial degree from list(range(degree))
    with bias (statistical distance from uniformity) which is
    O(2**-secpar)
    """
    if degree < 2:
        raise ValueError("Must have degree >= 2")
    return ceil(log2(degree)) + secpar


def bytes_for_index(secpar: int, degree: int) -> int:
    """
    Number of bytes to sample a monomial X**j with exponent j from
    list(range(degree)) with bias (statistical distance from uni-
    formity) which is O(2**-secpar)
    """
    return ceil(bits_for_index(secpar=secpar, degree=degree)/8)


def bits_for_fy_shuffle(secpar: int, degree: int) -> int:
    """
    Number of bits to Fisher-Yates shuffle list(range(degree)) with
    bias (statistical distance from uniformity) which is O(2**-secpar)
    """
    return degree * bits_for_index(secpar=secpar, degree=degree)


def bytes_for_fy_shuffle(secpar: int, degree: int) -> int:
    """
    Number of bytes to Fisher-Yates shuffle list(range(degree)) with
    bias (statistical distance from uniformity) which is O(2**-secpar)
    """
    return ceil(bits_for_fy_shuffle(secpar=secpar, degree=degree)/8)


def bits_for_bdd_poly(secpar: int, beta: int, degree: int, omega: int) -> int:
    """
    Number of bits to sample a polynomial which is a sum of omega
    monomials whose coefficients are in list(range(-beta,beta+1))
    with bias (statistical distance from uniformity) which is
    O(2**-secpar)
    """
    return omega * bits_for_bdd_coef(secpar=secpar, beta=beta) + bits_for_fy_shuffle(secpar=secpar, degree=degree)


def bytes_for_bdd_poly(secpar: int, beta: int, degree: int, omega: int) -> int:
    """
    Number of bytes to sample a polynomial which is a sum of omega
    monomials whose coefficients are in list(range(-beta,beta+1))
    with bias (statistical distance from uniformity) which is
    O(2**-secpar)
    """
    return ceil(bits_for_bdd_poly(secpar=secpar, beta=beta, degree=degree, omega=omega) / 8)


# Typing
SecretMat: type = Mat
PublicMat: type = Mat
SecretPoly: type = Poly
PublicPoly: type = Poly

