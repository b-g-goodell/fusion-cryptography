# algebra/polynomials.py
from typing import List, Tuple, Union
from api.ntt import ntt, cent, is_ntt_friendly, find_prou
from algebra.errors import (_MUST_BE_INT_ERR, _MUST_HAVE_PROU_ERR, _MUST_BE_LIST_ERR, _TYPE_MISMATCH_ERR,
                            _MODULUS_MISMATCH_ERR, _NTT_NOT_IMPLEMENTED_ERR, _NORM_NOT_IMPLEMENTED_ERR,
                            _MUL_BASE_NOT_IMPLEMENTED_ERR, _INVALID_REP_TYPE, _DEGREE_MISMATCH_ERR, _MUST_BE_TUPLE_ERR)


_POLYNOMIAL_REPRESENTATION_TYPES: List[str] = ['coefficient', 'ntt']
_POLYNOMIAL_REPRESENTATION_STR_PREFIX: str = '_PolynomialRepresentation'
_POLYNOMIAL_COEFFICIENT_REPRESENTATION_STR_PREFIX: str = '_PolynomialCoefficientRepresentation'
_POLYNOMIAL_NTT_REPRESENTATION_STR_PREFIX: str = '_PolynomialNTTRepresentation'


class PolynomialRepresentation:
    _mod: int
    vals: Tuple[int, ...]

    def __init__(self, mod: int, vals: Tuple[int, ...]):
        if not isinstance(vals, tuple):
            raise TypeError(_MUST_BE_TUPLE_ERR)
        elif not isinstance(mod, int) or not all(isinstance(i, int) for i in vals):
            raise TypeError(_MUST_BE_INT_ERR)
        elif not is_ntt_friendly(mod=mod, deg=len(vals)):
            raise ValueError(_MUST_HAVE_PROU_ERR)
        self._mod = mod
        self.vals = vals

    @property
    def mod(self) -> int:
        return self._mod

    @property
    def deg(self) -> int:
        return len(self.vals)

    @property
    def root_order(self) -> int:
        return 2*self.deg

    @property
    def root(self) -> int:
        return find_prou(mod=self.mod, deg=self.deg)

    @property
    def inv_root(self) -> int:
        return pow(base=self.root, exp=self.mod - 2, mod=self.mod)

    @property
    def coefs_norm_wght(self) -> Tuple[Tuple[int, ...], int, int]:
        raise NotImplementedError(_NORM_NOT_IMPLEMENTED_ERR)

    def __str__(self) -> str:
        return _POLYNOMIAL_REPRESENTATION_STR_PREFIX + f"(modulus={self.mod}, values={self.vals})"

    def __repr__(self) -> str:
        return self.__str__()

    def __eq__(self, other) -> bool:
        return (type(self) is type(other)
                and self.mod == other.mod
                and self.deg == other.deg
                and all((x - y) % self.mod == 0 for x, y in zip(self.vals, other.vals)))

    def __add__(self, other: 'PolynomialRepresentation') -> 'PolynomialRepresentation':
        if not type(self) is type(other):
            raise TypeError(_TYPE_MISMATCH_ERR)
        elif self.mod != other.mod:
            raise TypeError(_MODULUS_MISMATCH_ERR)
        elif self.deg != other.deg:
            raise TypeError(_DEGREE_MISMATCH_ERR)
        return PolynomialRepresentation(mod=self.mod, vals=tuple([cent(val=x + y, mod=self.mod) for x, y in zip(self.vals, other.vals)]))

    def __radd__(self, other: 'PolynomialRepresentation') -> 'PolynomialRepresentation':
        if other == 0:
            return self
        return self + other

    def __neg__(self) -> 'PolynomialRepresentation':
        return PolynomialRepresentation(mod=self.mod, vals=tuple([-x for x in self.vals]))

    def __sub__(self, other: 'PolynomialRepresentation') -> 'PolynomialRepresentation':
        return self + (-other)

    def __mul__(self, other: 'PolynomialRepresentation') -> 'PolynomialRepresentation':
        raise NotImplementedError(_MUL_BASE_NOT_IMPLEMENTED_ERR)

    def __rmul__(self, other: 'PolynomialRepresentation') -> 'PolynomialRepresentation':
        return self * other

    def transform(self) -> 'PolynomialRepresentation':
        raise NotImplementedError(_NTT_NOT_IMPLEMENTED_ERR)


def _is_mul_valid(left: PolynomialRepresentation, right: PolynomialRepresentation) -> bool:
    return type(left) is type(right) and left.mod == right.mod and left.deg == right.deg


class PolynomialCoefficientRepresentation(PolynomialRepresentation):
    def __init__(self, mod: int, vals: Tuple[int, ...]):
        super().__init__(mod=mod, vals=vals)

    @property
    def coefs_norm_wght(self) -> Tuple[Tuple[int, ...], int, int]:
        return self.vals, max(abs(x) for x in self.vals), sum(1 if x % self.mod != 0 else 0 for x in self.vals)

    def __str__(self) -> str:
        return _POLYNOMIAL_COEFFICIENT_REPRESENTATION_STR_PREFIX + f"(modulus={self.mod}, degree={self.deg}, root={self.root}, inv_root={self.inv_root}, values={self.vals})"

    def __add__(self, other: 'PolynomialCoefficientRepresentation') -> 'PolynomialCoefficientRepresentation':
        if not isinstance(other, PolynomialCoefficientRepresentation):
            raise TypeError(_TYPE_MISMATCH_ERR)
        return super().__add__(other=other)

    def __neg__(self) -> 'PolynomialCoefficientRepresentation':
        return PolynomialCoefficientRepresentation(mod=self.mod, vals=tuple([-x for x in self.vals]))

    def __mul__(self, other: 'PolynomialCoefficientRepresentation') -> 'PolynomialCoefficientRepresentation':
        if not _is_mul_valid(left=self, right=other):
            raise TypeError(_TYPE_MISMATCH_ERR)
        new_values: List[int] = [0 for _ in range(2 * self.deg)]
        for i, x in enumerate(self.vals):
            for j, y in enumerate(other.vals):
                new_values[i + j] += x * y
        new_values = [cent(val=x - y, mod=self.mod) for x, y in zip(new_values[:self.deg], new_values[self.deg:])]
        return PolynomialCoefficientRepresentation(mod=self.mod, vals=tuple(new_values))

    def transform(self) -> 'PolynomialNTTRepresentation':
        return PolynomialNTTRepresentation(mod=self.mod, vals=tuple(ntt(val=list(self.vals), mod=self.mod, inv_flag=False)))


class PolynomialNTTRepresentation(PolynomialRepresentation):
    def __init__(self, mod: int, vals: Tuple[int, ...]):
        super().__init__(mod=mod, vals=vals)

    @property
    def coefs_norm_wght(self) -> Tuple[Tuple[int, ...], int, int]:
        x: PolynomialCoefficientRepresentation = self.transform()
        return x.coefs_norm_wght

    def __str__(self) -> str:
        return _POLYNOMIAL_NTT_REPRESENTATION_STR_PREFIX + f"(modulus={self.mod}, degree={self.deg}, root={self.root}, inv_root={self.inv_root}, values={self.vals})"

    def __add__(self, other: 'PolynomialNTTRepresentation') -> 'PolynomialNTTRepresentation':
        if not isinstance(other, PolynomialNTTRepresentation):
            raise TypeError(_TYPE_MISMATCH_ERR)
        return super().__add__(other=other)

    def __neg__(self) -> 'PolynomialNTTRepresentation':
        return PolynomialNTTRepresentation(mod=self.mod, vals=tuple([-x for x in self.vals]))

    def __mul__(self, other: 'PolynomialNTTRepresentation') -> 'PolynomialNTTRepresentation':
        if other == 1:
            return self
        elif not _is_mul_valid(left=self, right=other):
            raise TypeError(_TYPE_MISMATCH_ERR)
        return PolynomialNTTRepresentation(mod=self.mod, vals=tuple([cent(val=x * y, mod=self.mod) for x, y in zip(self.vals, other.vals)]))

    def transform(self) -> PolynomialCoefficientRepresentation:
        return PolynomialCoefficientRepresentation(mod=self.mod, vals=tuple(ntt(val=list(self.vals), mod=self.mod, inv_flag=True)))


class _PolynomialFactory:
    @staticmethod
    def make(mod: int, vals: Tuple[int, ...], rep_flag: str) -> Union[PolynomialCoefficientRepresentation, PolynomialNTTRepresentation]:
        if rep_flag == _POLYNOMIAL_REPRESENTATION_TYPES[0]:
            return PolynomialCoefficientRepresentation(mod=mod, vals=vals)
        elif rep_flag == _POLYNOMIAL_REPRESENTATION_TYPES[1]:
            return PolynomialNTTRepresentation(mod=mod, vals=vals)
        raise ValueError(_INVALID_REP_TYPE + f": {rep_flag}")
