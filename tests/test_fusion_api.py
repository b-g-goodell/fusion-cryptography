from typing import List
from algebra.polynomials import transform
from fusion.fusion import PREFIX_PARAMETERS, GeneralMatrix, PolynomialNTTRepresentation, _hash_ch
from api.fusion_api import setup, generate_keys, generate_signature, aggregate_signatures, verify_aggregated_signature, Params, OneTimeKeyTuple, Signature, OneTimeVerificationKey, OneTimeSigningKey, SignatureChallenge
from secrets import randbelow

TEST_SECPAR = [128, 256]


def test_setup():
    for secpar in TEST_SECPAR:
        p: Params = setup(secpar=secpar)
        for k in p.__dict__:
            if k != 'public_challenge' and k != 'secpar':
                assert k in PREFIX_PARAMETERS[secpar]
                assert p.__dict__[k] == PREFIX_PARAMETERS[secpar][k]
        a: GeneralMatrix = p.public_challenge
        assert a.elem_class == PolynomialNTTRepresentation
        assert len(a.matrix) == PREFIX_PARAMETERS[secpar]["num_cols_sk"]
        for next_column in a.matrix:
            assert len(next_column) == PREFIX_PARAMETERS[secpar]["num_rows_sk"]
            for next_polynomial in next_column:
                assert isinstance(next_polynomial, PolynomialNTTRepresentation)
                assert next_polynomial.modulus == PREFIX_PARAMETERS[secpar]["modulus"]
                assert next_polynomial.degree == PREFIX_PARAMETERS[secpar]["degree"]
                assert next_polynomial.root == PREFIX_PARAMETERS[secpar]["root"]
                assert next_polynomial.inv_root == PREFIX_PARAMETERS[secpar]["inv_root"]
                assert next_polynomial.root_order == PREFIX_PARAMETERS[secpar]["root_order"]
                for next_coefficient in next_polynomial.values:
                    assert isinstance(next_coefficient, int)
                    assert -next_polynomial.halfmod <= next_coefficient <= next_polynomial.halfmod

def test_generate_keys():
    for secpar in TEST_SECPAR:
        # Setup
        params: Params = setup(secpar=secpar)

        num_keys: int = 2+randbelow(98)

        otks: List[OneTimeKeyTuple] = generate_keys(params=params, num_keys=num_keys)
        assert isinstance(otks, list)
        try:
            assert all(isinstance(otk, OneTimeKeyTuple) for otk in otks)
        except:
            for otk in otks:
                try:
                    assert isinstance(otk, OneTimeKeyTuple)
                except:
                    print("Non-otk found in output from key generation, " + str(type(otk)))
                    assert False
        for otk in otks:
            assert isinstance(otk, OneTimeKeyTuple)
            assert isinstance(otk.otsk, OneTimeSigningKey)
            assert isinstance(otk.otvk, OneTimeVerificationKey)

            # Unpack the key
            otsk: OneTimeSigningKey  # private
            otvk: OneTimeVerificationKey  # public
            otsk, otvk = otk.otsk, otk.otvk # unpack

            assert hasattr(otsk, "left_sk_hat")
            assert hasattr(otsk, "right_sk_hat")
            assert isinstance(otsk.left_sk_hat, GeneralMatrix)
            assert isinstance(otsk.right_sk_hat, GeneralMatrix)
            for each_sk_hat in [otsk.left_sk_hat, otsk.right_sk_hat]:
                assert len(each_sk_hat.matrix) == PREFIX_PARAMETERS[secpar]['num_rows_sk']
                for next_column in each_sk_hat.matrix:
                    assert len(next_column) == PREFIX_PARAMETERS[secpar]['num_cols_sk']
                    for next_polynomial in next_column:
                        assert isinstance(next_polynomial, PolynomialNTTRepresentation)
                        assert next_polynomial.degree == PREFIX_PARAMETERS[secpar]['degree']
                        assert next_polynomial.root == PREFIX_PARAMETERS[secpar]['root']
                        assert next_polynomial.inv_root == PREFIX_PARAMETERS[secpar]['inv_root']
                        assert next_polynomial.modulus == PREFIX_PARAMETERS[secpar]['modulus']
                        assert next_polynomial.root_order == PREFIX_PARAMETERS[secpar]['root_order']
                        for next_coefficient in next_polynomial.values:
                            assert isinstance(next_coefficient, int)
                            assert -next_polynomial.halfmod <= next_coefficient <= next_polynomial.halfmod

            inv_left_sk_hat: GeneralMatrix = GeneralMatrix(
                matrix=[[transform(x=f) for f in row] for row in otsk.left_sk_hat.matrix]
            )
            inv_right_sk_hat: GeneralMatrix = GeneralMatrix(
                matrix=[[transform(x=f) for f in row] for row in otsk.right_sk_hat.matrix]
            )
            assert inv_left_sk_hat.norm(p="infty") <= params.beta_sk
            assert inv_right_sk_hat.norm(p="infty") <= params.beta_sk
            assert inv_left_sk_hat.weight() <= params.omega_sk
            assert inv_right_sk_hat.weight() <= params.omega_sk

            assert hasattr(otvk, "left_vk_hat")
            assert hasattr(otvk, "right_vk_hat")
            assert isinstance(otvk.left_vk_hat, GeneralMatrix)
            assert isinstance(otvk.right_vk_hat, GeneralMatrix)

            assert params.public_challenge * otsk.left_sk_hat == otvk.left_vk_hat
            assert params.public_challenge * otsk.right_sk_hat == otvk.right_vk_hat


def test_generate_signature():
    for secpar in TEST_SECPAR:
        params: Params = setup(secpar=secpar)
        otk: OneTimeKeyTuple = generate_keys(params=params)
        otsk: OneTimeSigningKey
        otvk: OneTimeVerificationKey
        otsk, otvk = otk.otsk, otk.otvk
        message: str = "Hello, world!"
        ch: SignatureChallenge = _hash_ch(params=params, key=otvk, message=message)
        sig: Signature = generate_signature(params=params, key_pair=otk, message=message)
        assert isinstance(sig, Signature)
        assert isinstance(sig.signature_hat, GeneralMatrix)
        assert len(sig.signature_hat.matrix) == params.num_rows_sk
        assert len(sig.signature_hat.matrix[0]) == params.num_cols_sk
        assert all(
            isinstance(f, PolynomialNTTRepresentation) for row in sig.signature_hat.matrix for f in row
        )
        assert all(
            len(f.values) == params.degree
            for row in sig.signature_hat.matrix
            for f in row
        )
        target: GeneralMatrix = otvk.left_vk_hat * ch.c_hat + otvk.right_vk_hat
        observed: GeneralMatrix = params.public_challenge * sig.signature_hat
        assert target == observed

        inv_sig_hat: GeneralMatrix = GeneralMatrix(
            matrix=[[transform(f) for f in row] for row in sig.signature_hat.matrix]
        )
        sig_wt: int = inv_sig_hat.weight()
        assert sig_wt <= params.omega_vf_intermediate
        sig_norm: int = inv_sig_hat.norm(p="infty")
        assert sig_norm <= params.beta_vf_intermediate

# def test_aggregate_signatures():
#     for secpar in TEST_SECPAR:
#         params: Params = setup(secpar=secpar)
#         number_of_aggregands: int = 1+randbelow(params.capacity-1)
#         keys: List[OneTimeKeyTuple] = []
#         messages: List[str] = []
#         challenges: List[SignatureChallenge] = []
#         signatures: List[Signature] = []
#         for _ in range(number_of_aggregands):
#             keys += [generate_keys(params)]
#         keys: OneTimeKeyTuple = generate_keys(params=params)
#         sk: OneTimeSigningKey
#         vk: OneTimeVerificationKey
#         sk, vk = keys
#         message: str = "Hello, world!"
#         ch: SignatureChallenge = hash_ch(params=params, key=vk, message=message)
#         sig: Signature = generate_signature(params=params, key_pair=keys, message=message)
#         assert isinstance(sig, Signature)
#         assert isinstance(sig.signature_hat, GeneralMatrix)
#         assert len(sig.signature_hat.matrix) == params.num_rows_sk
#         assert len(sig.signature_hat.matrix[0]) == params.num_cols_sk
#         assert all(
#             isinstance(f, PolynomialNTTRepresentation) for row in sig.signature_hat.matrix for f in row
#         )
#         assert all(
#             len(f.values) == params.degree
#             for row in sig.signature_hat.matrix
#             for f in row
#         )
#         target: GeneralMatrix = vk.left_vk_hat * ch.c_hat + vk.right_vk_hat
#         observed: GeneralMatrix = params.public_challenge * sig.signature_hat
#         assert target == observed
#
#         inv_sig_hat: GeneralMatrix = GeneralMatrix(
#             matrix=[[transform(f) for f in row] for row in sig.signature_hat.matrix]
#         )
#         sig_wt: int = inv_sig_hat.weight()
#         assert sig_wt <= params.omega_vf_intermediate
#         sig_norm: int = inv_sig_hat.norm(p="infty")
#         assert sig_norm <= params.beta_vf_intermediate
#
# def test_verify_aggregated_signature():
#     pass