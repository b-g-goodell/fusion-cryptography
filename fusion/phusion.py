from typing import List, Dict, Any, Tuple, Iterator, Union
from hashlib import sha3_256, shake_256
from algebra.polynomials import Polynomial as Poly
from algebra.matrices import PolynomialMatrix as Mat, compute_lin_combo as _compute_lin_combo


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


class Params(object):
    secpar: int
    one: Poly
    pubchall: Mat
    bytes_to_sample_sk_coef: int
    bytes_to_sample_challenge: int
    bytes_to_sample_aggregators: int
    beta_vf: int
    omega_vf: int
    num_rows_sk: int
    num_cols_sk: int
    num_rows_pubchall: int
    num_cols_pubchall: int

    def __init__(self, secpar: int):
        self.secpar = secpar
        self.mod = IACR_SUGGESTED_PARAMS[secpar]['modulus']
        self.deg = IACR_SUGGESTED_PARAMS[secpar]['degree']
        self.one = Poly(
            mod=IACR_SUGGESTED_PARAMS[secpar]['modulus'],
            vals=tuple([1 for _ in range(IACR_SUGGESTED_PARAMS[secpar]['degree'])]),
            rep_flag='ntt')
        # self.pubchall = sample()


class OneTimeSigningKey(object):
    _val: Mat

    def __init__(self, vals: Mat):
        self._val = vals

    @property
    def val(self) -> Mat:
        return self._val


class OneTimeVerificationKey(object):
    _val: Mat

    def __init__(self, vals: Mat):
        self._val = vals

    @property
    def val(self) -> Mat:
        return self._val


class OneTimeKeys(object):
    _otsk: Tuple[OneTimeSigningKey, OneTimeSigningKey]
    _otvk: Tuple[OneTimeVerificationKey, OneTimeVerificationKey]

    def __init__(self, otsk: Tuple[OneTimeSigningKey, OneTimeSigningKey], otvk: Tuple[OneTimeVerificationKey, OneTimeVerificationKey]):
        self._otsk = otsk
        self._otvk = otvk

    @property
    def otsk(self) -> Tuple[OneTimeSigningKey, OneTimeSigningKey]:
        return self._otsk

    @property
    def otvk(self) -> Tuple[OneTimeVerificationKey, OneTimeVerificationKey]:
        return self._otvk

    @property
    def left_sk(self) -> Mat:
        return self._otsk[0].val

    @property
    def right_sk(self) -> Mat:
        return self._otsk[1].val

    @property
    def left_vk(self) -> Mat:
        return self._otvk[0].val

    @property
    def right_vk(self) -> Mat:
        return self._otvk[1].val


class Challenge(object):
    _val: Poly

    def __init__(self, val: Poly):
        self._val = val

    @property
    def val(self) -> Poly:
        return self._val


class Signature(object):
    _val: Mat

    def __init__(self, val: Mat):
        self._val = val

    @property
    def val(self) -> Mat:
        return self._val


class Aggregator(object):
    _val: Poly

    def __init__(self, val: Poly):
        self._val = val

    @property
    def val(self) -> Poly:
        return self._val


def fusion_setup(secpar: int) -> Params:
    return Params(**IACR_SUGGESTED_PARAMS[secpar])


def _bytes_to_sk(params: Params, randomness: bytes) -> Tuple[OneTimeSigningKey, OneTimeSigningKey]:
    # TODO: Fix this
    keys = []
    for i in range(2):
        keys += [[]]
        for j in range(params.num_rows_sk):
            keys[-1] += [[]]
            for k in range(params.num_cols_sk):
                keys[-1][-1] += [[]]
                for l in range(params.deg):
                    next_bytes, randomness = randomness[:params.bytes_to_sample_sk_coef], randomness[params.bytes_to_sample_sk_coef:]
                    keys[-1][-1][-1] += [int(next_bytes) % params.mod]
                next_poly = Poly(mod=params.mod, vals=keys[-1][-1][-1], rep_flag='coefficient')
    left_keys, right_keys = keys[0], keys[1]


def fusion_keygen(params: Params, randomness: bytes) -> OneTimeKeys:
    left_sk: OneTimeSigningKey = _bytes_to_sk(params=params, randomness=randomness[:len(randomness) // 2])
    left_vk: OneTimeVerificationKey = OneTimeVerificationKey(vals=params.pubchall * left_sk.val)
    rght_sk: OneTimeSigningKey = _bytes_to_sk(params=params, randomness=randomness[len(randomness) // 2:])
    rght_vk: OneTimeVerificationKey = OneTimeVerificationKey(vals=params.pubchall * rght_sk.val)
    return OneTimeKeys(otsk=(left_sk, rght_sk), otvk=(left_vk, rght_vk))


def _bytes_to_ch(params: Params, randomness: bytes) -> Challenge:
    # TODO: Finish
    pass


def _make_ch(params: Params, otvk: Tuple[OneTimeVerificationKey, OneTimeVerificationKey], message: bytes) -> Challenge:
    randomness: bytes = shake_256_wrapper(message=str(params)+str(otvk)+str(message), num_bytes=params.bytes_to_sample_challenge)
    return _bytes_to_ch(params=params, randomness=randomness)


def fusion_sign(params: Params, otk: OneTimeKeys, message: bytes) -> Signature:
    c: Challenge = _make_ch(params, otk.otvk, message)
    val: Mat = _compute_lin_combo(vectors=[x.val for x in otk.otsk], multipliers=[params.one, c.val])
    return Signature(val=val)


def _bytes_to_ags(params: Params, num_ags: int, randomness: bytes) -> List[Aggregator]:
    # TODO: Finish
    pass


def _sort_otvks_msgs(otvks: List[Tuple[OneTimeVerificationKey, OneTimeVerificationKey]], messages: List[bytes]) -> Tuple[List[Tuple[OneTimeVerificationKey, OneTimeVerificationKey]], List[bytes]]:
    return zip(*sorted(zip(otvks, messages), key=lambda x: x[0].val))


def _make_ch_and_ags(params: Params, otvks: List[Tuple[OneTimeVerificationKey, OneTimeVerificationKey]], messages: List[bytes]) -> Tuple[List[Challenge], List[Aggregator]]:
    n: int = len(otvks)
    srt_otvks, srt_msgs = _sort_otvks_msgs(otvks=otvks, messages=messages)
    srt_challs: List[Challenge] = []
    for otvk, message in otvks, messages:
        srt_challs += [_make_ch(params=params, otvk=otvk, message=message)]
    b: bytes = shake_256_wrapper(message=str(params)+str(srt_otvks)+str(srt_msgs)+str(srt_challs), num_bytes=params.bytes_to_sample_aggregators*n)
    aggregators: List[Aggregator] = _bytes_to_ags(params=params, num_ags=n, randomness=b)
    return srt_challs, aggregators


def _sort_otvks_msgs_sigs(otvks: List[Tuple[OneTimeVerificationKey, OneTimeVerificationKey]], messages: List[bytes], signatures: List[Signature]) -> Tuple[List[Tuple[OneTimeVerificationKey, OneTimeVerificationKey]], List[bytes], List[Signature]]:
    return zip(*sorted(zip(otvks, messages, signatures), key=lambda x: x[0].val))


def fusion_aggregate(params: Params, otvks: List[Tuple[OneTimeVerificationKey, OneTimeVerificationKey]], messages: List[bytes], signatures: List[Signature]) -> Signature:
    srt_otks, srt_msgs, srt_sigs = _sort_otvks_msgs_sigs(otvks=otvks, messages=messages, signatures=signatures)
    aggregators: List[Aggregator]
    _, aggregators = _make_ch_and_ags(params, srt_otks, srt_msgs)
    val: Mat = _compute_lin_combo(vectors=[x.val for x in srt_sigs], multipliers=[x.val for x in aggregators])
    return Signature(val=val)


def _make_target(otvks: List[Tuple[OneTimeVerificationKey, OneTimeVerificationKey]], challs: List[Challenge], aggregators: List[Aggregator]) -> Mat:
    vectors = [otvks[i][j].val for i in range(len(otvks)) for j in range(2)]
    multipliers = [aggregators[i].val if j == 0 else aggregators[i].val*challs[i].val for i in range(len(aggregators)) for j in range(2)]
    return _compute_lin_combo(vectors=vectors, multipliers=multipliers)


def fusion_verify(params: Params, otvks: List[Tuple[OneTimeVerificationKey, OneTimeVerificationKey]], messages: List[bytes], aggregate_signature: Signature) -> bool:
    srt_otvks, srt_msgs = _sort_otvks_msgs(otvks=otvks, messages=messages)
    challs: List[Challenge]
    aggregators: List[Aggregator]
    challs, aggregators = _make_ch_and_ags(params=params, otvks=srt_otvks, messages=srt_msgs)
    target: Mat = _make_target(otvks=otvks, challs=challs, aggregators=aggregators)
    evaluated_signature: Mat = params.pubchall * aggregate_signature.val
    c, n, w = aggregate_signature.val.coef_norm_wght
    return evaluated_signature == target and n <= params.beta_vf and w <= params.omega_vf