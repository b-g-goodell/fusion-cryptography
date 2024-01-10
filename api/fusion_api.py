from typing import List, Tuple
from fusion.fusion import fusion_setup, fusion_keygen, fusion_sign, aggregate, verify, Params, OneTimeKeyTuple, Signature, OneTimeVerificationKey, OneTimeSigningKey, SignatureChallenge

def setup(secpar: int) -> Params:
    """
    Setup for the fusion scheme with the given security parameter and random seed.
    """
    params = fusion_setup(secpar)
    return params

def generate_keys(params: Params, num_keys: int = 1) -> OneTimeKeyTuple | List[OneTimeKeyTuple]:
    """
    Generate a one-time signing key (SK) and one-time verification key (VK) pair.
    """
    if num_keys == 1:
        return fusion_keygen(params)
    return [fusion_keygen(params) for _ in range(num_keys)]

def generate_signature(params: Params, key_pair: OneTimeKeyTuple, message: str) -> Signature:
    """
    Generate a signature for a given message using the provided key pair.
    """
    signature = fusion_sign(params, key_pair, message)
    return signature

def aggregate_signatures(params: Params, keys: List[OneTimeVerificationKey], messages: List[str], signatures: List[Signature]) -> Signature:
    """
    Aggregate a list of signatures into a single signature using the provided keys and messages.
    """
    aggregate_signature = aggregate(params, keys, messages, signatures)
    return aggregate_signature

def verify_aggregated_signature(params: Params, keys: List[OneTimeVerificationKey], messages: List[str], aggregate_signature: Signature) -> Tuple[bool, str]:
    """
    Verify the validity of the aggregated signature given the keys and messages.
    """
    result, reason = verify(params, keys, messages, aggregate_signature)
    return result, reason