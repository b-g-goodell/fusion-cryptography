from crypto.fusion.errors import _ParameterMismatchError
from crypto.fusion.params import Params
from algebra.polynomials import _Polynomial as _Poly
from algebra.matrices import _PolynomialMatrix as _PolyMat
from algebra.errors import _DimensionMismatchError, _TypeMismatchError


class SecretMatrix(_PolyMat):
    vals: _PolyMat

    def __init__(self, vals: _PolyMat):
        if not isinstance(vals, _PolyMat):
            raise _DimensionMismatchError
        super().__init__(vals=vals.vals)


class PublicMatrix(_PolyMat):
    vals: _PolyMat

    def __init__(self, vals: _PolyMat):
        if not isinstance(vals, _PolyMat):
            raise _DimensionMismatchError
        super().__init__(vals=vals.vals)


class SecretPolynomial(_Poly):
    val: _Poly

    def __init__(self, val: _Poly):
        if not isinstance(val, _Poly):
            raise _DimensionMismatchError
        super().__init__(mod=val.modulus, vals=val.vals, rep_flag='ntt')


class PublicPolynomial(_Poly):
    val: _Poly

    def __init__(self, val: _Poly):
        if not isinstance(val, _Poly):
            raise _DimensionMismatchError
        super().__init__(mod=val.modulus, vals=val.vals, rep_flag='ntt')


class OneTimeSigningKey:
    left: SecretMatrix
    rght: SecretMatrix

    def __init__(self, left: SecretMatrix, rght: SecretMatrix):
        if not isinstance(left, SecretMatrix) or not isinstance(rght, SecretMatrix):
            raise _ParameterMismatchError
        self.left = left
        self.rght = rght

    def __eq__(self, other: 'OneTimeSigningKey') -> bool:
        return self.left == other.left and self.rght == other.rght

    def __str__(self) -> str:
        return self.__class__.__name__ + f"(left={self.left}, right={self.rght})"

    def __repr__(self) -> str:
        return self.__str__()


class OneTimeVerificationKey:
    left: PublicMatrix
    rght: PublicMatrix

    def __init__(self, left: PublicMatrix, rght: PublicMatrix):
        if not isinstance(left, PublicMatrix) or not isinstance(rght, PublicMatrix):
            raise _ParameterMismatchError
        self.left = left
        self.rght = rght

    def __eq__(self, other: 'OneTimeVerificationKey') -> bool:
        return self.left == other.left and self.rght == other.rght

    def __str__(self) -> str:
        return self.__class__.__name__ + f"(left={self.left}, right={self.rght})"

    def __repr__(self) -> str:
        return self.__str__()


class OneTimeKeyPair:
    otsk: OneTimeSigningKey
    otvk: OneTimeVerificationKey

    def __init__(self, otsk: OneTimeSigningKey, otvk: OneTimeVerificationKey):
        if not isinstance(otsk, OneTimeSigningKey) or not isinstance(otvk, OneTimeVerificationKey):
            raise _TypeMismatchError
        self.otsk = otsk
        self.otvk = otvk

    def __eq__(self, other: 'OneTimeKeyPair') -> bool:
        return self.otsk == other.otsk and self.otvk == other.otvk

    def __str__(self) -> str:
        return self.__class__.__name__ + f"(otsk={self.otsk}, otvk={self.otvk})"

    def __repr__(self) -> str:
        return self.__str__()


class SignatureChallenge:
    val: PublicPolynomial

    def __init__(self, val: PublicPolynomial):
        if not isinstance(val, PublicPolynomial):
            raise _TypeMismatchError
        self.val = val

    def __eq__(self, other: 'SignatureChallenge') -> bool:
        return self.val == other.val

    def __str__(self) -> str:
        return self.__class__.__name__ + f"(val={self.val})"

    def __repr__(self) -> str:
        return self.__str__()

    @classmethod
    def make_challenge(cls, otvk: OneTimeVerificationKey, msg: bytes) -> 'SignatureChallenge':
