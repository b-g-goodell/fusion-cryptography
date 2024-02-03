from typing import Dict, Union, Tuple, List
from algebra.polynomials import Polynomial as Poly
from algebra.matrices import PolynomialMatrix as Mat
from fusion.sampling import sample_matrix_by_ntt

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


# Typing
SecretMat: type = Mat
PublicMat: type = Mat
SecretPoly: type = Poly
PublicPoly: type = Poly
OneTimeSigningKey: type = Tuple[SecretMat, ...]
OneTimeVerificationKey: type = Tuple[PublicMat, ...]
OneTimeKeyTuple: type = Tuple[OneTimeSigningKey, OneTimeVerificationKey]


def parse_bytes_as_public_challenge(deg: int, modulus: int, rows: int, cols: int, randomness: bytes) -> PublicMat:
    values: List[List[Poly]] = []
    for i in range(rows):
        values += [[]]
        for j in range(cols):
            values[-1] += [parse_bytes_as_poly]


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
    public_challenge: PublicMat
    sign_pre_hash_dst: str
    sign_hash_dst: str
    agg_xof_dst: str
    bytes_for_one_coef_bdd_by_beta_ch: int
    bytes_for_one_coef_bdd_by_beta_ag: int
    bytes_for_fy_shuffle: int
    bytes_per_sig_chall_poly: int
    bytes_per_agg_coef_poly: int

    def __init__(self, secpar: int, seed: bytes):
        if secpar not in IACR_SUGGESTED_PARAMS:
            raise ValueError("Invalid security parameter.")
        self.secpar = secpar
        self.seed = seed
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
        self.public_challenge = parse_bytes_as_public_challenge(IACR_SUGGESTED_PARAMS[secpar]["public_challenge"])
        # self.public_challenge = sample_matrix_by_ntt(modulus=self.modulus, degree=self.degree, num_rows=self.num_rows_pub_challenge, num_cols=self.num_rows_sk)

    def __str__(self) -> str:
        return f"Params(secpar={self.secpar}, capacity={self.capacity}, modulus={self.modulus}, degree={self.degree}, root_order={self.root_order}, root={self.root}, inv_root={self.inv_root}, num_rows_pub_challenge={self.num_rows_pub_challenge}, num_rows_sk={self.num_rows_sk}, num_cols_sk={self.num_cols_sk}, beta_sk={self.beta_sk}, beta_ch={self.beta_ch}, beta_ag={self.beta_ag}, beta_vf={self.beta_vf}, omega_sk={self.omega_sk}, omega_ch={self.omega_ch}, omega_ag={self.omega_ag}, omega_vf={self.omega_vf}, public_challenge={str(self.public_challenge)}, sign_pre_hash_dst={self.sign_pre_hash_dst}, sign_hash_dst={self.sign_hash_dst}, agg_xof_dst={self.agg_xof_dst}, bytes_for_one_coef_bdd_by_beta_ch={self.bytes_for_one_coef_bdd_by_beta_ch}, bytes_for_one_coef_bdd_by_beta_ag={self.bytes_for_one_coef_bdd_by_beta_ag}, bytes_for_fy_shuffle={self.bytes_for_fy_shuffle})"

    def __eq__(self, other):
        return self.__dict__ == other.__dict__
