from hashlib import shake_256, sha3_256
from math import ceil, log2
from typing import List, Dict, Tuple, Union
from sampling.sampling import sample_matrix_by_coefs, sample_matrix_by_ntt
from crypto.fusion.errors import (_MUST_BE_MATRIX_ERR, _ELEM_CLASS_MISMATCH_ERR, _NORM_TOO_LARGE_ERR, _KEYS_NOT_VALID_ERR,
                                  _WGHT_TOO_LARGE_ERR, _MUST_BE_POLY_ERR, _PARAMS_MISMATCH_ERR, _CHALL_NOT_VALID_ERR,
                                  _LENGTH_MISMATCH, _AGG_COEFS_NOT_VALID_ERR, _MUST_BE_PARAMS_ERR)
from algebra.polynomials import _Polynomial as Poly
from algebra.matrices import _PolynomialMatrix as Matrix


# Some global constants
PRIME: int = 2147465729
DEGREES: Dict[int, int] = {128: 2**6, 256: 2**8}
RANKS: Dict[int, int] = {128: 195, 256: 83}
CAPACITIES: Dict[int, int] = {128: 1796, 256: 2818}
ROOTS: Dict[int, int] = {128: 23584283, 256: 3337519}
SIGN_PRE_HASH_DSTS: Dict[int, str] = {128: '0000', 256: '0001'}
SIGN_HASH_DSTS: Dict[int, str] = {128: '0010', 256: '0011'}
AGG_XOF_DSTS: Dict[int, str] = {128: '0100', 256: '0101'}
CH_WTS: Dict[int, int] = {128: 27, 256: 60}
AG_WTS: Dict[int, int] = {128: 35, 256: 60}
CH_BDS: Dict[int, int] = {128: 3, 256: 1}
AG_BDS: Dict[int, int] = {128: 2, 256: 1}
IACR_SUGGESTED_PARAMS: Dict[int, Dict[str, Union[int, str]]] = {secpar: {
    "capacity": CAPACITIES[secpar], "modulus": PRIME, "degree": DEGREES[secpar],
    "root_order": 2 * DEGREES[secpar], "root": ROOTS[secpar], "inv_root": pow(ROOTS[secpar], PRIME - 2, PRIME),
    "num_rows_pub_challenge": 1, "num_rows_sk": RANKS[secpar], "num_rows_vk": 1,
    "num_cols_sk": 1, "num_cols_vk": 1, "sign_pre_hash_dst": SIGN_PRE_HASH_DSTS[secpar],
    "sign_hash_dst": SIGN_HASH_DSTS[secpar], "agg_xof_dst": AGG_XOF_DSTS[secpar], "beta_sk": 52, "beta_ch": 1,
    "beta_ag": 1, "omega_sk": DEGREES[secpar], "omega_ch": CH_WTS[secpar], "omega_ag": AG_WTS[secpar]
} for secpar in [128, 256]}

# __str__ prefixes
OTSK_STR_PREFIX: str = "OneTimeSigningKey"
OTVK_STR_PREFIX: str = "OneTimeVerificationKey"
OTK_STR_PREFIX: str = "OneTimeKeyTuple"
SIG_CHALL_PREFIX: str = "SignatureChallenge"
MSG_PREFIX: str = "Message"
SIG_PREFIX: str = "Signature"
AGG_COEF_PREFIX: str = "AggregationCoefficient"
AGG_SIG_PREFIX: str = "AggregateSignature"
KMC_PREFIX = "MessageOTVKChallengeTuple"
MOTVKCSignatureTuple_PREFIX = "MessageOTVKChallengeSignatureTuple"
AGG_COEFS_MATRIX_PREFIX: str = "AggCoefsMatrix"
AGGREGATION_PACKAGE_PREFIX: str = "AggregationPackage"

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


# simple functions for counting bits and bytes
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
class Params(object):
    secpar: int
    capacity: int
    modulus: int
    degree: int
    root_order: int
    root: int
    inv_root: int
    num_rows_pub_challenge: int
    num_rows_sk: int
    num_cols_sk: int
    beta_sk: int
    beta_ch: int
    beta_ag: int
    beta_vf_intermediate: int
    beta_vf: int
    omega_sk: int
    omega_ch: int
    omega_ag: int
    omega_vf_intermediate: int
    omega_vf: int
    public_challenge: Matrix
    sign_pre_hash_dst: str
    sign_hash_dst: str
    agg_xof_dst: str
    bytes_for_one_coef_bdd_by_beta_ch: int
    bytes_for_one_coef_bdd_by_beta_ag: int
    bytes_for_fy_shuffle: int
    bytes_per_sig_chall_poly: int
    bytes_per_agg_coef_poly: int

    def __init__(self, secpar: int):
        if secpar not in IACR_SUGGESTED_PARAMS:
            raise ValueError("Invalid security parameter.")
        self.secpar = secpar
        self.capacity = IACR_SUGGESTED_PARAMS[secpar]["capacity"]
        self.modulus = IACR_SUGGESTED_PARAMS[secpar]["modulus"]
        self.degree = IACR_SUGGESTED_PARAMS[secpar]["degree"]
        self.root_order = IACR_SUGGESTED_PARAMS[secpar]["root_order"]
        self.root = IACR_SUGGESTED_PARAMS[secpar]["root"]
        self.inv_root = IACR_SUGGESTED_PARAMS[secpar]["inv_root"]
        self.num_rows_pub_challenge = IACR_SUGGESTED_PARAMS[secpar]["num_rows_pub_challenge"]
        self.num_rows_sk = IACR_SUGGESTED_PARAMS[secpar]["num_rows_sk"]
        # self.num_rows_vk = PREFIX_PARAMETERS[secpar]["num_rows_vk"] == self.num_rows_pub_challenge
        self.num_cols_sk = IACR_SUGGESTED_PARAMS[secpar]["num_cols_sk"]
        # self.num_cols_vk = PREFIX_PARAMETERS[secpar]["num_cols_vk"] == self.num_cols_sk
        self.sign_pre_hash_dst = IACR_SUGGESTED_PARAMS[secpar]["sign_pre_hash_dst"]
        self.sign_hash_dst = IACR_SUGGESTED_PARAMS[secpar]["sign_hash_dst"]
        self.agg_xof_dst = IACR_SUGGESTED_PARAMS[secpar]["agg_xof_dst"]
        self.beta_sk = IACR_SUGGESTED_PARAMS[secpar]["beta_sk"]
        self.beta_ch = IACR_SUGGESTED_PARAMS[secpar]["beta_ch"]
        self.beta_ag = IACR_SUGGESTED_PARAMS[secpar]["beta_ag"]
        self.beta_vf_intermediate = IACR_SUGGESTED_PARAMS[secpar]["beta_vf_intermediate"]
        self.beta_vf = IACR_SUGGESTED_PARAMS[secpar]["beta_vf"]
        self.omega_sk = IACR_SUGGESTED_PARAMS[secpar]["omega_sk"]
        self.omega_ch = IACR_SUGGESTED_PARAMS[secpar]["omega_ch"]
        self.omega_ag = IACR_SUGGESTED_PARAMS[secpar]["omega_ag"]
        self.omega_vf_intermediate = IACR_SUGGESTED_PARAMS[secpar]["omega_vf_intermediate"]
        self.omega_vf = IACR_SUGGESTED_PARAMS[secpar]["omega_vf"]
        self.public_challenge = sample_matrix_by_ntt(modulus=self.modulus, degree=self.degree, num_rows=self.num_rows_pub_challenge, num_cols=self.num_rows_sk)

    def __str__(self) -> str:
        return f"Params(secpar={self.secpar}, capacity={self.capacity}, modulus={self.modulus}, degree={self.degree}, root_order={self.root_order}, root={self.root}, inv_root={self.inv_root}, num_rows_pub_challenge={self.num_rows_pub_challenge}, num_rows_sk={self.num_rows_sk}, num_cols_sk={self.num_cols_sk}, beta_sk={self.beta_sk}, beta_ch={self.beta_ch}, beta_ag={self.beta_ag}, beta_vf={self.beta_vf}, omega_sk={self.omega_sk}, omega_ch={self.omega_ch}, omega_ag={self.omega_ag}, omega_vf={self.omega_vf}, public_challenge={str(self.public_challenge)}, sign_pre_hash_dst={self.sign_pre_hash_dst}, sign_hash_dst={self.sign_hash_dst}, agg_xof_dst={self.agg_xof_dst}, bytes_for_one_coef_bdd_by_beta_ch={self.bytes_for_one_coef_bdd_by_beta_ch}, bytes_for_one_coef_bdd_by_beta_ag={self.bytes_for_one_coef_bdd_by_beta_ag}, bytes_for_fy_shuffle={self.bytes_for_fy_shuffle})"

    def __eq__(self, other):
        return self.__dict__ == other.__dict__


class OneTimeSigningKey(object):
    _params: Params
    _left: Matrix
    _rght: Matrix

    def __init__(self, params: Params, left: Matrix, rght: Matrix):
        if not isinstance(left, Matrix) or not isinstance(rght, Matrix):
            raise TypeError(_MUST_BE_MATRIX_ERR)
        elif not left.elem_class is rght.elem_class or not rght.elem_class is Poly:
            raise ValueError(_ELEM_CLASS_MISMATCH_ERR)
        elif params.public_challenge.rows != params.num_rows_pub_challenge or params.public_challenge.cols != params.num_rows_sk:
            raise ValueError(_DIMENSION_MISMATCH_ERR)
        elif not left.elem_class is rght.elem_class or not rght.elem_class is Poly:
            raise ValueError(_ELEM_CLASS_MISMATCH_ERR)
        elif left.rows != rght.rows or left.rows != params.num_rows_sk or left.cols != rght.cols or left.cols != params.num_cols_sk:
            raise ValueError(_DIMENSION_MISMATCH_ERR)
        elif any(f.deg != params.degree for next_row in left.vals + rght.vals for f in next_row):
            raise TypeError(_DEGREE_MISMATCH_ERR)
        for next_row in left.vals + rght.vals:
            for f in next_row:
                if not isinstance(f, Poly):
                    raise TypeError(_MUST_BE_POLY_ERR)
                elif f._ntt_rep.modulus != params.modulus:
                    raise ValueError(_MODULUS_MISMATCH_ERR)
                elif f._ntt_rep.deg != params.degree:
                    raise ValueError(_DEGREE_MISMATCH_ERR)
                elif f._ntt_rep.root_order != params.root_order:
                    raise ValueError(_ROOT_ORDER_MISMATCH_ERR)
                elif f._ntt_rep.root != params.root:
                    raise ValueError(_ROOT_MISMATCH_ERR)
                elif f._ntt_rep.inv_root != params.inv_root:
                    raise ValueError(_INV_ROOT_MISMATCH_ERR)
                _, n, w = f.coef_norm_wght
                if n > params.beta_sk:
                    raise ValueError(_NORM_TOO_LARGE_ERR)
                elif w > params.omega_sk:
                    raise ValueError(_WGHT_TOO_LARGE_ERR)
        self._params = params
        self._left = left
        self._rght = rght

    def __str__(self):
        # Do not print your key.
        return OTSK_STR_PREFIX + f"(params={self.params},left={str(self.left)}, rght={str(self.rght)})"

    @property
    def params(self) -> Params:
        return self._params

    @property
    def left(self) -> Matrix:
        return self._left

    @property
    def rght(self) -> Matrix:
        return self._rght

    def __eq__(self, other):
        return (self.params == other.params and
                (self.left - other.left) % self.params.modulus == 0 and
                (self.rght - other.rght) % self.params.modulus == 0)


class OneTimeVerificationKey(object):
    _params: Params
    _left: Matrix
    _rght: Matrix

    def __init__(self, params: Params, left: Matrix, rght: Matrix):
        if not isinstance(params, Params):
            raise TypeError(_MUST_BE_PARAMS_ERR)
        elif not isinstance(left, Matrix) or not isinstance(rght, Matrix):
            raise TypeError(_MUST_BE_MATRIX_ERR)
        elif params.public_challenge.rows != params.num_rows_pub_challenge or params.public_challenge.cols != params.num_rows_sk:
            raise ValueError(_DIMENSION_MISMATCH_ERR)
        elif not left.elem_class is rght.elem_class or not rght.elem_class is Poly:
            raise ValueError(_ELEM_CLASS_MISMATCH_ERR)
        elif left.rows != rght.rows or left.rows != params.num_rows_pub_challenge or left.cols != rght.cols or left.cols != params.num_cols_sk:
            raise ValueError(DIMENSION_MISMATCH_ERR)
        for next_row in left.vals + rght.vals:
            for next_poly in next_row:
                if not isinstance(next_poly, Poly):
                    raise TypeError(_MUST_BE_POLY_ERR)
                elif next_poly._ntt_rep.modulus != params.modulus:
                    raise ValueError(_MODULUS_MISMATCH_ERR)
                elif next_poly._ntt_rep.deg != params.degree or len(next_poly._ntt_rep._vals) != params.degree:
                    raise ValueError(_DEGREE_MISMATCH_ERR)
                elif next_poly._ntt_rep.root_order != params.root_order:
                    raise ValueError(_ROOT_ORDER_MISMATCH_ERR)
                elif next_poly._ntt_rep.root != params.root:
                    raise ValueError(_ROOT_MISMATCH_ERR)
                elif next_poly._ntt_rep.inv_root != params.inv_root:
                    raise ValueError(_INV_ROOT_MISMATCH_ERR)
        self._params = params
        self._left = left
        self._rght = rght

    def __str__(self):
        # Do not print your key.
        return OTSK_STR_PREFIX + f"(params={self.params},left={str(self.left)}, rght={str(self.rght)})"

    @property
    def params(self) -> Params:
        return self._params

    @property
    def left(self) -> Matrix:
        return self._left

    @property
    def rght(self) -> Matrix:
        return self._rght

    def __eq__(self, other):
        return (self.params == other.params and
                (self.left - other.left) % self.params.modulus == 0 and
                (self.rght - other.rght) % self.params.modulus == 0)


class OneTimeKey(object):
    _otsk: OneTimeSigningKey
    _otvk: OneTimeVerificationKey

    def __init__(self, params: Params, otsk: OneTimeSigningKey, otvk: OneTimeVerificationKey):
        if otsk.params != otvk.params:
            raise ValueError(_PARAMS_MISMATCH_ERR)
        elif ((params.public_challenge * otsk.left - otvk.left) % params.modulus
              or (params.public_challenge * otsk.rght - otvk.rght) % params.modulus):
            raise ValueError(_KEYS_NOT_VALID_ERR)
        self._otsk = otsk
        self._otvk = otvk

    @property
    def params(self) -> Params:
        return self.otvk.params

    @property
    def otsk(self) -> OneTimeSigningKey:
        return self._otsk

    @property
    def otvk(self) -> OneTimeVerificationKey:
        return self._otvk

    def __str__(self):
        return OTK_STR_PREFIX+f"(params={self.params},otsk={self.otsk},otvk={self.otvk})"

    def __eq__(self, other):
        return self.otsk == other.otsk and self.otvk == other.otvk


class Challenge(object):
    _params: Params
    _val: Poly

    def __init__(self, params: Params, val: Poly):
        if not isinstance(val, Poly):
            raise TypeError(_MUST_BE_POLY_ERR)
        elif val._ntt_rep.modulus != params.modulus:
            raise ValueError(_MODULUS_MISMATCH_ERR)
        elif val._ntt_rep.deg != params.degree or len(val._ntt_rep._vals) != params.degree:
            raise ValueError(_DEGREE_MISMATCH_ERR)
        elif val._ntt_rep.root_order != params.root_order:
            raise ValueError(_ROOT_ORDER_MISMATCH_ERR)
        elif val._ntt_rep.root != params.root:
            raise ValueError(_ROOT_MISMATCH_ERR)
        elif val._ntt_rep.inv_root != params.inv_root:
            raise ValueError(_INV_ROOT_MISMATCH_ERR)
        elif not isinstance(params, Params):
            raise TypeError(_MUST_BE_PARAMS_ERR)
        c, n, w = val.coef_norm_wght
        if n > params.beta_ch:
            raise ValueError(_NORM_TOO_LARGE_ERR)
        elif w > params.omega_ch:
            raise ValueError(_WGHT_TOO_LARGE_ERR)
        self._params = params
        self._val = val

    def __str__(self):
        return SIG_CHALL_PREFIX+f"(params={self.params},val={str(self.val)})"

    def __eq__(self, other):
        return self.chall == other.chall

    @property
    def params(self) -> Params:
        return self._params

    @property
    def val(self) -> Poly:
        return self._val


# Supporting internal functions
def _pre_hash_msg_to_int_digest(params: Params, msg: str) -> int:
    # Apply sha3-256 for hash-then-sign
    dst_and_msg: bytes = str.encode(
        params.sign_pre_hash_dst
        + ","
        + msg
    )
    pre_hashed_message: bytes = sha3_256_wrapper(dst_and_msg)
    pre_hashed_message_as_int: int = int.from_bytes(pre_hashed_message, byteorder="little")
    return pre_hashed_message_as_int


def _hash_vk_and_pre_hashed_msg_to_bytes_digest(params: Params, otvk: OneTimeVerificationKey, pre_hashed_msg: int, num_bytes: int) -> bytes:
    dst_and_otk_and_prehashed_msg: bytes = str.encode(
        params.sign_hash_dst
        + ","
        + str(otvk)
        + ","
        + str(pre_hashed_msg)
    )
    return shake_256_wrapper(dst_and_otk_and_prehashed_msg, num_bytes)


def _bytes_digest_to_polys(params: Params, beta: int, omega: int, b: bytes, num_polys: int) -> List[Poly]:
    if omega < 1:
        raise ValueError("Cannot decode bytes to no coefficients.")
    elif beta < 1:
        raise ValueError("Cannot decode bytes to non-positive coefficients.")
    elif num_polys < 1:
        raise ValueError("Cannot decode to no polynomials.")
    coef_bytes: int = bytes_for_bdd_coef(secpar=params.secpar, beta=beta)
    idx_bytes: int = bytes_for_index(secpar=params.secpar, degree=params.degree)
    bytes_per_bdd_poly: int = bytes_for_bdd_poly(secpar=params.secpar, beta=beta, degree=params.degree, omega=omega)
    if len(b) < num_polys * bytes_per_bdd_poly:
        raise ValueError(
            f"Too few bytes to decode polynomial. Expected {bytes_per_bdd_poly} but got {len(b)}"
        )
    result: List[Poly] = []
    next_b: bytes
    remaining_b: bytes = b
    while remaining_b or next_b:
        next_b, remaining_b = (
            remaining_b[:bytes_per_bdd_poly],
            remaining_b[bytes_per_bdd_poly:]
        )
        # Generate the coefficients
        coefficients: List[int] = []
        next_bytes: bytes
        remaining_bytes: bytes = remaining_b
        for i in range(omega):
            next_bytes, remaining_bytes = (
                remaining_bytes[:coef_bytes],
                remaining_bytes[coef_bytes:],
            )
            coefficients += [
                (int.from_bytes(next_bytes, byteorder="little") % (2*beta + 1)) - beta
            ]
        while len(coefficients) < params.degree:
            coefficients += [0]

        # Fisher-Yates shuffle if necessary
        # This part can be optimized a lot
        if omega < params.degree:
            for i in range(params.degree - 1, omega, -1):
                # convert next_few_bytes to an integer modulo i+1
                next_bytes, remaining_bytes = (
                    remaining_bytes[:idx_bytes],
                    remaining_bytes[idx_bytes:],
                )
                j: int = int.from_bytes(next_bytes, byteorder="little") % (i + 1)
                coefficients[i], coefficients[j] = coefficients[j], coefficients[i]
        result += [Poly(modulus=params.modulus, values=coefficients, representation='coefficient')]
    return result


def _make_sig_challenge(params: Params, otvk: OneTimeVerificationKey, msg: str) -> Challenge:
    num_bytes_needed: int = bytes_for_bdd_poly(secpar=params.secpar, beta=params.beta_ch, degree=params.degree, omega=params.omega_ch)
    prehashed_msg: int = _pre_hash_msg_to_int_digest(params=params, msg=msg)
    hashed_key_and_prehashed_msg: bytes = _hash_vk_and_pre_hashed_msg_to_bytes_digest(params=params, otvk=otvk, pre_hashed_msg=prehashed_msg, num_bytes=num_bytes_needed)
    c: Poly = _bytes_digest_to_polys(params=params, beta=params.beta_ch, omega=params.omega_ch, b=hashed_key_and_prehashed_msg, num_polys=1)[0]
    return Challenge(params=params, val=c)


class KeyMessageChallenge:
    _otvk: OneTimeVerificationKey
    _msg: str
    _chall: Challenge

    def __init__(self, msg: str, otvk: OneTimeVerificationKey, chall: Challenge):
        if chall != _make_sig_challenge(params=otvk.params, otvk=otvk, msg=msg):
            raise ValueError(_CHALL_NOT_VALID_ERR)
        self._otvk = otvk
        self._msg = msg
        self._chall = chall

    def __str__(self):
        return KMC_PREFIX+ f"(params={self.params},chall={str(self.chall)})"

    def __eq__(self, other):
        return self.msg == other.msg and self.otvk == other.otvk and self.chall == other.chall

    @property
    def params(self) -> Params:
        return self.otvk.params

    @property
    def otvk(self) -> OneTimeVerificationKey:
        return self._otvk

    @property
    def msg(self) -> str:
        return self._msg

    @property
    def chall(self) -> Challenge:
        return self._chall


class Signature(object):
    _signature: Matrix

    def __init__(self, signature: Matrix):
        for next_row in signature.vals:
            for next_poly in next_row:
                _, n, w = next_poly.coef_norm_wght
                if n > signature.params.beta_vf_intermediate:
                    raise ValueError(_NORM_TOO_LARGE_ERR)
                elif w > signature.params.omega_vf_intermediate:
                    raise ValueError(_WGHT_TOO_LARGE_ERR)
        self.params = params
        self._signature = signature

    def __str__(self):
        return SIG_PREFIX+f"(params={self.params},signature={str(self._signature)})"

    def __eq__(self, other):
        return self._signature == other._signature and self.params == other.params

    @property
    def params(self) -> Params:
        return self._params


class SignedMessage(object):
    kmc: KeyMessageChallenge
    sig: Signature

    def __init__(self, params: Params, motvkc: KeyMessageChallenge, sig: Signature):
        if params != motvkc.params or params != sig.params:
            raise ValueError(_PARAMS_MISMATCH_ERR)
        self.params = params
        self.kmc = motvkc
        self.sig = sig

    def __str__(self):
        return MOTVKCSignatureTuple_PREFIX + f"(params={self.params},motvkc={str(self.kmc)},sig={str(self.sig)})"

    def __eq__(self, other):
        return self.params == other.params and self.kmc == other.kmc and self.sig == other.sig

    @property
    def params(self) -> Params:
        return self.params


class AggregationCoefficient(object):
    _params: Params
    alpha: Poly

    def __init__(self, params: Params, alpha: Poly):
        c, n, w = alpha.coef_norm_wght
        if n > params.beta_ag:
            raise ValueError(_NORM_TOO_LARGE_ERR)
        elif w > params.omega_ag:
            raise ValueError(_WGHT_TOO_LARGE_ERR)
        self.params = params
        self.alpha = alpha

    def __str__(self):
        return AGG_COEF_PREFIX+f"(params={self.params},alpha={self.alpha})"

    def __eq__(self, other):
        return self.alpha == other.alpha and self.params == other.params

    @property
    def params(self) -> Params:
        return self._params


class AggregationMatrix(object):
    # Todo: finish type checking
    agg_coefs: Matrix

    def __init__(self, params: Params, agg_coefs: Matrix):
        elem_class: type = type(None)
        for each_agg_coef in agg_coefs:
            if elem_class is type(None):
                elem_class = type(each_agg_coef)
            elif not isinstance(each_agg_coef, elem_class):
                raise TypeError(_ELEM_CLASS_MISMATCH_ERR)
            c, n, w = each_agg_coef.coef_norm_wght
            if n > params.beta_ag:
                raise ValueError(_NORM_TOO_LARGE_ERR)
            elif w > params.omega_ag:
                raise ValueError(_WGHT_TOO_LARGE_ERR)
        self.params = params
        self.agg_coefs = agg_coefs

    def __str__(self):
        return AGG_COEFS_MATRIX_PREFIX + f"(params={self.params},agg_coefs={self.agg_coefs})"

    def __eq__(self, other):
        return self.params == other.params and self.agg_coefs == other.agg_coefs

    @property
    def params(self) -> Params:
        return self.agg_coefs[0][0].params


def _make_aggregation_coefficients(params: Params, srt_otvks: List[OneTimeVerificationKey],
                                   srt_pre_hashed_msgs: List[int], srt_challs: List[Challenge]) -> List[AggregationCoefficient]:
    num_keys: int = len(srt_otvks)
    num_bytes_needed_per_poly: int = bytes_for_bdd_poly((params.secpar, params.beta_ag, params.degree, params.omega_ag))
    num_bytes_needed: int = num_keys * num_bytes_needed_per_poly
    bytes_digest: bytes = _hash_vks_and_pre_hashed_msgs_and_challs_to_bytes_digest(params=params, otvks=srt_otvks, pre_hashed_msgs=srt_pre_hashed_msgs, challs=srt_challs, num_bytes=num_bytes_needed)
    alphas: List[Poly] = _bytes_digest_to_polys(params=params, beta=params.beta_ag, degree=params.degree, omega=params.omega_ag, b=bytes_digest, num_polys=num_keys)
    return [AggregationCoefficient(alpha=alpha) for alpha in alphas]


class AggregationPackage(object):
    params: Params
    motvkcs: List[KeyMessageChallenge]
    sigs: List[Signature]
    agg_coefs: AggregationMatrix

    def __init__(self, params: Params, motvkcs: List[KeyMessageChallenge], sigs: List[Signature], agg_coefs: AggregationMatrix):
        if not all(params == x.params for x in motvkcs) or not all(params == x.params for x in sigs) or params != agg_coefs.params:
            raise ValueError(_PARAMS_MISMATCH_ERR)
        elif len(motvkcs) != len(sigs) or len(motvkcs) != len(agg_coefs.agg_coefs):
            raise ValueError(_LENGTH_MISMATCH)
        srt_stuff = sorted([x for x in motvkcs], key=lambda x: x.otvk)
        if agg_coefs != _make_aggregation_coefficients(params=params, srt_otvks=[x.otvk for x in srt_stuff], srt_pre_hashed_msgs=[x.msg for x in srt_stuff], srt_challs=[x.chall for x in srt_stuff]):
            raise ValueError(_AGG_COEFS_NOT_VALID_ERR)
        self.params = params
        self.motvkcs = motvkcs
        self.sigs = sigs
        self.agg_coefs = agg_coefs

    def __str__(self):
        return AGGREGATION_PACKAGE_PREFIX + f"(params={self.params},motvkcs={self.motvkcs},sigs={self.sigs},agg_coefs={self.agg_coefs})"

    def __eq__(self, other):
        return self.params == other.params and self.motvkcs == other.motvkcs and self.sigs == other.sigs and self.agg_coefs == other.agg_coefs


class AggregateSignature(object):
    aggregate_signature: Matrix

    def __init__(self, params: Params, aggregate_signature: Matrix):
        if not isinstance(aggregate_signature, Matrix):
            raise TypeError(_MUST_BE_MATRIX_ERR)
        c, n, w = aggregate_signature.coefs_norm_weight
        if n > params.beta_vf:
            raise ValueError(_NORM_TOO_LARGE_ERR)
        elif w > params.omega_vf:
            raise ValueError(_WGHT_TOO_LARGE_ERR)
        self.aggregate_signature = aggregate_signature

    def __str__(self):
        return AGG_SIG_PREFIX+f"(aggregate_signature={str(self.aggregate_signature)})"


def _verify_signature(params: Params, otvk: OneTimeVerificationKey, msg: str, sig: Signature) -> Tuple[bool, str]:
    c, n, w = sig._signature.coefs_norm_weight
    if n > params.beta_vf_intermediate:
        return False, "Norm of purported signature too large"
    elif w > params.omega_vf_intermediate:
        return False, "Hamming weight of purported signature too large"
    c: Challenge = _make_sig_challenge(params=params, otvk=otvk, msg=msg)
    target_image: Matrix = otvk.left_vk_hat + otvk.right_vk_hat * c
    if params.public_challenge * sig != target_image:
        return False, "Signature image does not match target image."
    return True, "Signature valid."


def _hash_vks_and_pre_hashed_msgs_and_challs_to_bytes_digest(params: Params, otvks: List[OneTimeVerificationKey], pre_hashed_msgs: List[int], challs: List[Challenge], num_bytes: int) -> bytes:
    dst_and_vks_and_pre_hashed_msgs_and_challs: bytes = str.encode(
        params.agg_xof_dst
        + ","
        + str(list(zip(otvks, pre_hashed_msgs, challs)))
    )
    return shake_256_wrapper(dst_and_vks_and_pre_hashed_msgs_and_challs, num_bytes)


# Fusion algorithms
def fusion_setup(secpar: int) -> Params:
    return Params(secpar=secpar)


def fusion_keygen(params: Params) -> OneTimeKey:
    left_key_coefs: Matrix = sample_matrix_by_coefs(mod=params.modulus, deg=params.degree,
                                                    num_rows=params.num_rows_sk, num_cols=params.num_cols_sk,
                                                    norm=params.beta_sk, wght=params.omega_sk)
    right_key_coefs: Matrix = sample_matrix_by_coefs(mod=params.modulus, deg=params.degree,
                                                     num_rows=params.num_rows_sk, num_cols=params.num_cols_sk,
                                                     norm=params.beta_sk, wght=params.omega_sk)
    left_sk_hat: Matrix = Matrix(
        matrix=[[transform(y) for y in z] for z in left_key_coefs.vals]
    )
    right_sk_hat: Matrix = Matrix(
        matrix=[[transform(y) for y in z] for z in right_key_coefs.vals]
    )
    otsk: OneTimeSigningKey = OneTimeSigningKey(left_sk_hat=left_sk_hat, right_sk_hat=right_sk_hat)
    left_vk_hat: Matrix = params.public_challenge * left_sk_hat
    right_vk_hat: Matrix = params.public_challenge * right_sk_hat
    otvk: OneTimeVerificationKey = OneTimeVerificationKey(left_vk_hat=left_vk_hat, right_vk_hat=right_vk_hat)
    return OneTimeKey(otsk=otsk, otvk=otvk)


def fusion_sign(params: Params, otk: OneTimeKey, msg: str) -> Signature:
    c_hat: Challenge = _make_sig_challenge(params=params, otvk=otk.otvk, msg=msg)
    return Signature(signature_hat=otk.otsk.left_sk * c_hat.c_hat + otk.otsk.rght_sk)


# def fusion_aggregate(params: Params, otvks: List[OneTimeVerificationKey], msgs: List[str], sigs: List[Signature]) -> AggregateSignature:
#     pre_hashed_messages: List[int] = [_pre_hash_msg_to_int_digest(params=params, msg=msg) for msg in msgs]
#     challs: List[SignatureChallenge] = [_make_sig_challenge(params=params, otvk=otvk, msg=msg) for otvk, msg in zip(otvks, msgs)]
#     sorted_input: List[Tuple[OneTimeVerificationKey, str, Signature, int, SignatureChallenge]] = sorted(
#         list(zip(otvks, msgs, sigs, pre_hashed_messages, challs)), key=lambda x: str(x[0])
#     )
#     sorted_otvks: List[Signature] = [x[0] for x in sorted_input]
#     sorted_signatures: List[Signature] = [x[2] for x in sorted_input]
#     sorted_prehashed_messages: List[int] = [x[3] for x in sorted_input]
#     sorted_challs: List[SignatureChallenge] = [x[4] for x in sorted_input]
#
#     aggregation_coefficients: List[AggregationCoefficient] = _make_aggregation_coefficients(params=params,
#                                                                                             srt_otvks=sorted_otvks,
#                                                                                             srt_pre_hashed_msgs=sorted_prehashed_messages,
#                                                                                             srt_challs=sorted_challs)
#     aggregate_signature_hat_values: Matrix = (sorted_signatures[0].signature_hat * aggregation_coefficients[0].alpha_hat)
#     for next_alpha, next_sig in zip(aggregation_coefficients[1:], sorted_signatures[1:]):
#         aggregate_signature_hat_values += next_sig.signature_hat * next_alpha.alpha_hat
#     return Signature(signature_hat=aggregate_signature_hat_values)

#
# def fusion_verify(
#     params: Params,
#     keys: List[OneTimeVerificationKey],
#     messages: List[str],
#     aggregate_signature: Signature,
# ) -> Tuple[bool, str]:
#     if len(keys) > params.capacity:
#         return False, f"Too many keys."
#     elif len(keys) != len(messages):
#         return False, f"Number of keys and messages must be equal."
#     coef_rep_agg_sig: Matrix = Matrix(
#         matrix=[[transform(z) for z in y] for y in aggregate_signature.signature_hat]
#     )
#     sorted_input = sorted(zip(keys, messages), key=lambda x: str(x[0]))
#     sorted_vks: List[OneTimeVerificationKey] = [x[0] for x in sorted_input]
#     sorted_challs: List[SignatureChallenge] = [
#         _hash_ch(params=params, key=next_vk, message=next_m)
#         for next_vk, next_m in sorted_input
#     ]
#     aggregation_coefficients: List[AggregationCoefficient] = _make_agg_coefs(params=params,
#                                                                              keys=[x[0] for x in sorted_input],
#                                                                              messages=[x[1] for x in sorted_input])
#     tmp: Matrix = sorted_vks[0].left_vk_hat * sorted_challs[0].c_hat
#     tmp += sorted_vks[0].right_vk_hat
#     target: Matrix = (
#         sorted_vks[0].left_vk_hat * sorted_challs[0].c_hat + sorted_vks[0].right_vk_hat
#     ) * aggregation_coefficients[0].alpha_hat
#     for next_alpha, next_vk, next_chall in zip(
#         aggregation_coefficients[1:], sorted_vks[1:], sorted_challs[1:]
#     ):
#         target += (
#             next_vk.left_vk_hat * next_chall.c_hat + next_vk.right_vk_hat
#         ) * next_alpha.alpha_hat
#     observed: Matrix = (
#         params.public_challenge * aggregate_signature.signature_hat
#     )
#     for a, b in zip(target.matrix, observed.matrix):
#         for c, d in zip(a, b):
#             if c != d:
#                 return False, f"Target doesn't match image of aggregate signature."
#     if any(
#         z.norm(p="infty") > params.beta_vf for y in coef_rep_agg_sig.matrix for z in y
#     ):
#         return False, f"Norm of aggregate signature too large."
#     elif any(z.weight() > params.omega_vf for y in coef_rep_agg_sig.matrix for z in y):
#         return False, f"Weight of aggregate signature too large."
#     return True, ""
