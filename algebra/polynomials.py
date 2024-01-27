from typing import Dict as _Dict, List as _List, Tuple as _Tuple
from api.ntt import ntt, cent, has_prou, is_prou, is_root_inverse
from algebra.errors import (_MUST_BE_INT_ERR, _MUST_HAVE_PROU_ERR, _MUST_BE_CORRECT_ROOT_ERR,
                            _MUST_BE_CORRECT_INVERSE_ROOT_ERR, _MUST_BE_LIST_ERR, _DEGREE_MISMATCH_ERR,
                            _TYPE_MISMATCH_ERR, _MODULUS_MISMATCH_ERR, _ROOT_MISMATCH_ERR, _NTT_NOT_IMPLEMENTED_ERR,
                            _NORM_NOT_IMPLEMENTED_ERR, _MUL_BASE_NOT_IMPLEMENTED_ERR)

_CACHED_ROOT_ORDERS: _Dict[int, int] = {}
_CACHED_HALFMODS: _Dict[int, int] = {}
_CACHED_LOGMODS: _Dict[int, int] = {}


class _PolynomialRepresentation:
    modulus: int
    degree: int
    root: int
    inv_root: int
    values: _List[int]

    def __init__(self, modulus: int, degree: int, root: int, inv_root: int, values: _List[int]):
        if not isinstance(modulus, int) or not isinstance(degree, int) or not isinstance(root, int) or not isinstance(inv_root, int) or not all(isinstance(i, int) for i in values):
            raise TypeError(_MUST_BE_INT_ERR)
        elif not isinstance(values, list):
            raise TypeError(_MUST_BE_LIST_ERR)
        elif len(values) != degree:
            raise TypeError(_DEGREE_MISMATCH_ERR)
        elif not has_prou(mod=modulus, deg=degree):
            raise ValueError(_MUST_HAVE_PROU_ERR)
        elif not is_prou(root=root, mod=modulus, deg=degree):
            raise ValueError(_MUST_BE_CORRECT_ROOT_ERR)
        elif not is_root_inverse(root=root, inv_root=inv_root, mod=modulus):
            raise ValueError(_MUST_BE_CORRECT_INVERSE_ROOT_ERR)
        self.modulus = modulus
        self.degree = degree
        self.root = root
        self.inv_root = inv_root
        self.values = values

    def __str__(self):
        return f"PolynomialRepresentation(modulus={self.modulus}, degree={self.degree}, root={self.root}, inv_root={self.inv_root}, values={self.values})"

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        if not isinstance(other, _PolynomialRepresentation):
            return False
        elif self.modulus != other.modulus:
            return False
        elif self.degree != other.degree:
            return False
        elif self.root != other.root:
            return False
        elif self.inv_root != other.inv_root:
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
        elif self.degree != other.degree:
            raise TypeError(_DEGREE_MISMATCH_ERR)
        elif self.root != other.root:
            raise TypeError(_ROOT_MISMATCH_ERR)
        return _PolynomialCoefficientRepresentation(
            modulus=self.modulus,
            degree=self.degree,
            root=self.root,
            inv_root=self.inv_root,
            values=[
                cent(val=x + y, mod=self.modulus)
                for x, y in zip(self.values, other.values)])

    def __radd__(self, other):
        return self + other

    def __neg__(self):
        return _PolynomialRepresentation(
            modulus=self.modulus,
            degree=self.degree,
            root=self.root,
            inv_root=self.inv_root,
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

    @property
    def halfmod(self) -> int:
        if self.modulus not in _CACHED_HALFMODS:
            _CACHED_HALFMODS[self.modulus] = self.modulus // 2
        return _CACHED_HALFMODS[self.modulus]

    @property
    def logmod(self) -> int:
        if self.modulus not in _CACHED_LOGMODS:
            _CACHED_LOGMODS[self.modulus] = self.modulus.bit_length()
        return _CACHED_LOGMODS[self.modulus]

    @property
    def root_order(self) -> int:
        if self.degree not in _CACHED_ROOT_ORDERS:
            _CACHED_ROOT_ORDERS[self.degree] = 2 * self.degree
        return _CACHED_ROOT_ORDERS[self.degree]


def _is_mul_valid(left: _PolynomialRepresentation, right: _PolynomialRepresentation):
    if not isinstance(left, type(right)) or not isinstance(right, type(left)):
        raise TypeError(_TYPE_MISMATCH_ERR)
    elif left.modulus != right.modulus:
        raise TypeError(_MODULUS_MISMATCH_ERR)
    elif left.degree != right.degree or len(left.values) != len(right.values) or len(left.values) != left.degree:
        raise TypeError(_DEGREE_MISMATCH_ERR)
    elif left.root != right.root:
        raise TypeError(_ROOT_MISMATCH_ERR)


class _PolynomialCoefficientRepresentation(_PolynomialRepresentation):
    def __init__(
        self,
        modulus: int,
        degree: int,
        root: int,
        inv_root: int,
        values: _List[int]
    ):
        super().__init__(
            modulus=modulus,
            degree=degree,
            root=root,
            inv_root=inv_root,
            values=values)

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
            degree=self.degree,
            root=self.root,
            inv_root=self.inv_root,
            values=[-(x % self.modulus) for x in self.values])

    def __mul__(self, other):
        if other == 0:
            return 0
        elif other == 1:
            return self
        _is_mul_valid(left=self, right=other)
        values: _List[int] = [0 for _ in range(2 * self.degree)]
        for i, x in enumerate(self.values):
            for j, y in enumerate(other.values):
                values[i + j] += x * y
        values = [cent(val=x - y, mod=self.modulus) for x, y in zip(values[:self.degree], values[self.degree:])]
        return _PolynomialCoefficientRepresentation(
            modulus=self.modulus,
            degree=self.degree,
            root=self.root,
            inv_root=self.inv_root,
            values=values)

    def __rmul__(self, other):
        return self * other

    def transform(self) -> '_PolynomialNTTRepresentation':
        return _PolynomialNTTRepresentation(
            modulus=self.modulus,
            degree=self.degree,
            root=self.root,
            inv_root=self.inv_root,
            values=ntt(val=self.values, mod=self.modulus, inv_flag=False))

    @property
    def norm_weight(self) -> _Tuple[int, int]:
        return max(abs(x) for x in self.values), sum(1 if x % self.modulus != 0 else 0 for x in self.values)


class _PolynomialNTTRepresentation(_PolynomialRepresentation):
    def __init__(
        self,
        modulus: int,
        degree: int,
        root: int,
        inv_root: int,
        values: _List[int]
    ):
        super().__init__(
            modulus=modulus,
            degree=degree,
            root=root,
            inv_root=inv_root,
            values=values)

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
            degree=self.degree,
            root=self.root,
            inv_root=self.inv_root,
            values=[-(x % self.modulus) for x in self.values])

    def __mul__(self, other):
        if other == 0:
            return 0
        elif other == 1:
            return self
        _is_mul_valid(left=self, right=other)
        return _PolynomialNTTRepresentation(
            modulus=self.modulus,
            degree=self.degree,
            root=self.root,
            inv_root=self.inv_root,
            values=[cent(val=x * y, mod=self.modulus) for x, y in zip(self.values, other.values)])

    def transform(self) -> _PolynomialCoefficientRepresentation:
        return _PolynomialCoefficientRepresentation(
            modulus=self.modulus,
            degree=self.degree,
            root=self.root,
            inv_root=self.inv_root,
            values=ntt(val=self.values, mod=self.modulus, inv_flag=True))

    @property
    def norm_weight(self) -> _Tuple[int, int]:
        x: _PolynomialCoefficientRepresentation = self.transform()
        return x.norm_weight


class _PolynomialFactory:
    @staticmethod
    def create_representation(modulus: int, degree: int, root: int, inv_root: int, values: _List[int], representation_type: str):
        if representation_type == 'coefficient':
            return _PolynomialCoefficientRepresentation(modulus=modulus, degree=degree, root=root, inv_root=inv_root, values=values)
        elif representation_type == 'ntt':
            # Additional logic can go here to handle NTT-specific data
            return _PolynomialNTTRepresentation(modulus=modulus, degree=degree, root=root, inv_root=inv_root, values=values)
        else:
            raise ValueError(f"Unknown representation type: {representation_type}")
