from typing import List as _List, Tuple as _Tuple
from api.ntt import ntt, cent, is_ntt_friendly, find_prou
from algebra.errors import (_MUST_BE_INT_ERR, _MUST_HAVE_PROU_ERR, _MUST_BE_LIST_ERR, _TYPE_MISMATCH_ERR,
                            _MODULUS_MISMATCH_ERR, _NTT_NOT_IMPLEMENTED_ERR, _NORM_NOT_IMPLEMENTED_ERR,
                            _MUL_BASE_NOT_IMPLEMENTED_ERR)


class _PolynomialRepresentation:
    _modulus: int
    _values: _List[int]

    def __init__(self, modulus: int, values: _List[int]):
        if not isinstance(values, list):
            raise TypeError(_MUST_BE_LIST_ERR)
        elif not isinstance(modulus, int) or not all(isinstance(i, int) for i in values):
            raise TypeError(_MUST_BE_INT_ERR)
        elif not is_ntt_friendly(mod=modulus, deg=len(values)):
            raise ValueError(_MUST_HAVE_PROU_ERR)
        self._modulus = modulus
        self._values = values

    @property
    def modulus(self) -> int:
        return self._modulus

    @property
    def values(self) -> _List[int]:
        return self._values

    @property
    def degree(self) -> int:
        return len(self.values)

    @property
    def root_order(self) -> int:
        return 2*self.degree

    @property
    def root(self) -> int:
        return find_prou(mod=self.modulus, deg=self.degree)

    @property
    def inv_root(self) -> int:
        return pow(base=self.root, exp=self.modulus-2, mod=self.modulus)

    def __str__(self):
        return f"PolynomialRepresentation(modulus={self.modulus}, values={self.values})"

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        if not isinstance(other, _PolynomialRepresentation):
            return False
        elif self.modulus != other.modulus:
            return False
        return all(
            (x - y) % self.modulus == 0
            for x, y in zip(self.values, other.values))

    def __add__(self, other):
        if other == 0:
            return self
        elif not isinstance(other, _PolynomialRepresentation):
            raise TypeError(_TYPE_MISMATCH_ERR)
        elif self.modulus != other.modulus:
            raise TypeError(_MODULUS_MISMATCH_ERR)
        return _PolynomialRepresentation(
            modulus=self.modulus,
            values=[
                cent(val=x + y, mod=self.modulus)
                for x, y in zip(self.values, other.values)])

    def __radd__(self, other):
        return self + other

    def __neg__(self):
        return _PolynomialRepresentation(
            modulus=self.modulus,
            values=[-(x % self.modulus) for x in self.values])

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        raise NotImplementedError(_MUL_BASE_NOT_IMPLEMENTED_ERR)

    def __rmul__(self, other):
        return self * other

    def transform(self):
        raise NotImplementedError(_NTT_NOT_IMPLEMENTED_ERR)

    @property
    def norm_weight(self) -> _Tuple[int, int]:
        raise NotImplementedError(_NORM_NOT_IMPLEMENTED_ERR)


def _is_mul_valid(left: _PolynomialRepresentation, right: _PolynomialRepresentation)  -> bool:
    return isinstance(left, type(right)) and isinstance(right, type(left)) and left.modulus == right.modulus and len(left.values) == len(right.values)


class _PolynomialCoefficientRepresentation(_PolynomialRepresentation):
    def __init__(self,modulus: int, values: _List[int]):
        super().__init__(modulus=modulus, values=values)

    def __str__(self):
        return f"PolynomialCoefficientRepresentation(modulus={self.modulus}, degree={self.degree}, root={self.root}, inv_root={self.inv_root}, values={self.values})"

    def __eq__(self, other):
        if not isinstance(other, _PolynomialCoefficientRepresentation):
            return False
        return super().__eq__(other=other)

    def __add__(self, other):
        if other == 0:
            return self
        elif not isinstance(other, _PolynomialCoefficientRepresentation):
            raise TypeError(_TYPE_MISMATCH_ERR)
        return super().__add__(other=other)

    def __neg__(self):
        return _PolynomialCoefficientRepresentation(
            modulus=self.modulus,
            values=[-(x % self.modulus) for x in self.values])

    def __mul__(self, other):
        if other == 0:
            return 0
        elif other == 1:
            return self
        elif not _is_mul_valid(left=self, right=other):
            raise TypeError(_TYPE_MISMATCH_ERR)
        new_values: _List[int] = [0 for _ in range(2 * self.degree)]
        for i, x in enumerate(self.values):
            for j, y in enumerate(other.values):
                new_values[i + j] += x * y
        new_values = [cent(val=x - y, mod=self.modulus) for x, y in zip(new_values[:self.degree], new_values[self.degree:])]
        return _PolynomialCoefficientRepresentation(modulus=self.modulus, values=new_values)

    def __rmul__(self, other):
        return self * other

    def transform(self) -> '_PolynomialNTTRepresentation':
        return _PolynomialNTTRepresentation(
            modulus=self.modulus,
            values=ntt(val=self.values, mod=self.modulus, inv_flag=False))

    @property
    def norm_weight(self) -> _Tuple[int, int]:
        return max(abs(x) for x in self.values), sum(1 if x % self.modulus != 0 else 0 for x in self.values)


class _PolynomialNTTRepresentation(_PolynomialRepresentation):
    def __init__(self, modulus: int, values: _List[int]):
        super().__init__(modulus=modulus, values=values)

    def __str__(self):
        return f"PolynomialNTTRepresentation(modulus={self.modulus}, degree={self.degree}, root={self.root}, inv_root={self.inv_root}, values={self.values})"

    def __eq__(self, other):
        if not isinstance(other, _PolynomialNTTRepresentation):
            return False
        return super().__eq__(other=other)

    def __add__(self, other):
        if other == 0:
            return self
        elif not isinstance(other, _PolynomialNTTRepresentation):
            raise TypeError(_TYPE_MISMATCH_ERR)
        return super().__add__(other=other)

    def __neg__(self):
        return _PolynomialNTTRepresentation(
            modulus=self.modulus,
            values=[-(x % self.modulus) for x in self.values])

    def __mul__(self, other):
        if other == 0:
            return 0
        elif other == 1:
            return self
        _is_mul_valid(left=self, right=other)
        return _PolynomialNTTRepresentation(
            modulus=self.modulus,
            values=[cent(val=x * y, mod=self.modulus) for x, y in zip(self.values, other.values)])

    def transform(self) -> _PolynomialCoefficientRepresentation:
        return _PolynomialCoefficientRepresentation(
            modulus=self.modulus,
            values=ntt(val=self.values, mod=self.modulus, inv_flag=True))

    @property
    def norm_weight(self) -> _Tuple[int, int]:
        x: _PolynomialCoefficientRepresentation = self.transform()
        return x.norm_weight


class _PolynomialFactory:
    @staticmethod
    def create_representation(modulus: int, values: _List[int], representation_type: str):
        if representation_type == 'coefficient':
            return _PolynomialCoefficientRepresentation(modulus=modulus, values=values)
        elif representation_type == 'ntt':
            return _PolynomialNTTRepresentation(modulus=modulus, values=values)
        raise ValueError(f"Unknown representation type: {representation_type}")
