# algebra/polynomials.py
from typing import List, Tuple, Union
from api.ntt import ntt, cent, is_ntt_friendly, find_prou
from algebra.errors import (_MUST_BE_INT_ERR, _MUST_HAVE_PROU_ERR, _MUST_BE_LIST_ERR, _TYPE_MISMATCH_ERR,
                            _MODULUS_MISMATCH_ERR, _NTT_NOT_IMPLEMENTED_ERR, _NORM_NOT_IMPLEMENTED_ERR,
                            _MUL_BASE_NOT_IMPLEMENTED_ERR, _INVALID_REP_TYPE_ERR, _DEGREE_MISMATCH_ERR, _MUST_BE_TUPLE_ERR)


_POLYNOMIAL_REPRESENTATION_TYPES: List[str] = ['coefficient', 'ntt']
_POLYNOMIAL_REPRESENTATION_STR_PREFIX: str = '_PolynomialRepresentation'
_POLYNOMIAL_COEFFICIENT_REPRESENTATION_STR_PREFIX: str = '_PolynomialCoefficientRepresentation'
_POLYNOMIAL_NTT_REPRESENTATION_STR_PREFIX: str = '_PolynomialNTTRepresentation'
_POLYNOMIAL_REPRESENTATION_STR_PREFIX: str = 'Polynomial'


class _PolynomialRepresentation:
    _mod: int
    _vals: Tuple[int, ...]

    def __init__(self, mod: int, vals: Tuple[int, ...]):
        if not isinstance(vals, tuple):
            raise TypeError(_MUST_BE_TUPLE_ERR)
        elif not isinstance(mod, int) or not all(isinstance(i, int) for i in vals):
            raise TypeError(_MUST_BE_INT_ERR)
        elif not is_ntt_friendly(mod=mod, deg=len(vals)):
            raise ValueError(_MUST_HAVE_PROU_ERR)
        self._mod = mod
        self._vals = vals

    @property
    def mod(self) -> int:
        return self._mod

    @property
    def deg(self) -> int:
        return len(self._vals)

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

    @property
    def vals(self) -> Tuple[int, ...]:
        return self._vals

    def __str__(self) -> str:
        return _POLYNOMIAL_REPRESENTATION_STR_PREFIX + f"(modulus={self.mod}, values={self._vals})"

    def __eq__(self, other) -> bool:
        return (type(self) is type(other)
                and self.mod == other.mod
                and self.deg == other.deg
                and all((x - y) % self.mod == 0 for x, y in zip(self._vals, other._val)))

    def __add__(self, other: '_PolynomialRepresentation') -> '_PolynomialRepresentation':
        if not type(self) is type(other):
            raise TypeError(_TYPE_MISMATCH_ERR)
        elif self.mod != other.mod:
            raise TypeError(_MODULUS_MISMATCH_ERR)
        elif self.deg != other.deg:
            raise TypeError(_DEGREE_MISMATCH_ERR)
        return _PolynomialRepresentation(mod=self.mod, vals=tuple([cent(val=x + y, mod=self.mod) for x, y in zip(self._vals, other._vals)]))

    def __radd__(self, other: '_PolynomialRepresentation') -> '_PolynomialRepresentation':
        if other == 0:
            return self
        return self + other

    def __neg__(self) -> '_PolynomialRepresentation':
        return _PolynomialRepresentation(mod=self.mod, vals=tuple([-x for x in self._vals]))

    def __sub__(self, other: '_PolynomialRepresentation') -> '_PolynomialRepresentation':
        return self + (-other)

    def __mul__(self, other: '_PolynomialRepresentation') -> '_PolynomialRepresentation':
        raise NotImplementedError(_MUL_BASE_NOT_IMPLEMENTED_ERR)

    def __rmul__(self, other: '_PolynomialRepresentation') -> '_PolynomialRepresentation':
        return self * other

    def __mod__(self, other: int) -> '_PolynomialRepresentation':
        if not isinstance(other, int):
            raise TypeError(_MUST_BE_INT_ERR)
        return _PolynomialRepresentation(mod=other, vals=tuple([cent(val=x, mod=other) for x in self._vals]))

    def transform(self) -> '_PolynomialRepresentation':
        raise NotImplementedError(_NTT_NOT_IMPLEMENTED_ERR)

    @property
    def coefs_norm_wght(self) -> Tuple[Tuple[int, ...], int, int]:
        raise NotImplementedError(_NTT_NOT_IMPLEMENTED_ERR)


def _is_mul_valid(left: _PolynomialRepresentation, right: _PolynomialRepresentation) -> bool:
    return type(left) is type(right) and left.mod == right.mod and left.deg == right.deg


class _PolynomialCoefficientRepresentation(_PolynomialRepresentation):
    def __init__(self, mod: int, vals: Tuple[int, ...]):
        super().__init__(mod=mod, vals=vals)

    def __str__(self) -> str:
        return _POLYNOMIAL_COEFFICIENT_REPRESENTATION_STR_PREFIX + f"(modulus={self.mod}, degree={self.deg}, root={self.root}, inv_root={self.inv_root}, values={self._vals})"

    def __add__(self, other: '_PolynomialCoefficientRepresentation') -> '_PolynomialCoefficientRepresentation':
        if not isinstance(other, _PolynomialCoefficientRepresentation):
            raise TypeError(_TYPE_MISMATCH_ERR)
        return super().__add__(other=other)

    def __neg__(self) -> '_PolynomialCoefficientRepresentation':
        return _PolynomialCoefficientRepresentation(mod=self.mod, vals=tuple([-x for x in self._vals]))

    def __mul__(self, other: '_PolynomialCoefficientRepresentation') -> '_PolynomialCoefficientRepresentation':
        if not _is_mul_valid(left=self, right=other):
            raise TypeError(_TYPE_MISMATCH_ERR)
        new_values: List[int] = [0 for _ in range(2 * self.deg)]
        for i, x in enumerate(self._vals):
            for j, y in enumerate(other._vals):
                new_values[i + j] += x * y
        new_values = [cent(val=x - y, mod=self.mod) for x, y in zip(new_values[:self.deg], new_values[self.deg:])]
        return _PolynomialCoefficientRepresentation(mod=self.mod, vals=tuple(new_values))

    def __mod__(self, other: int) -> '_PolynomialCoefficientRepresentation':
        if not isinstance(other, int):
            raise TypeError(_MUST_BE_INT_ERR)
        return _PolynomialCoefficientRepresentation(mod=other, vals=tuple([cent(val=x, mod=other) for x in self._vals]))

    def transform(self) -> '_PolynomialNTTRepresentation':
        return _PolynomialNTTRepresentation(mod=self.mod, vals=tuple(ntt(val=list(self._vals), mod=self.mod, inv_flag=False)))

    @property
    def coefs_norm_wght(self) -> Tuple[Tuple[int, ...], int, int]:
        return self._vals, max(abs(x) for x in self._vals), sum(1 if x % self.mod != 0 else 0 for x in self._vals)


class _PolynomialNTTRepresentation(_PolynomialRepresentation):
    def __init__(self, mod: int, vals: Tuple[int, ...]):
        super().__init__(mod=mod, vals=vals)

    def __str__(self) -> str:
        return _POLYNOMIAL_NTT_REPRESENTATION_STR_PREFIX + f"(modulus={self.mod}, degree={self.deg}, root={self.root}, inv_root={self.inv_root}, values={self._vals})"

    def __add__(self, other: '_PolynomialNTTRepresentation') -> '_PolynomialNTTRepresentation':
        if not isinstance(other, _PolynomialNTTRepresentation):
            raise TypeError(_TYPE_MISMATCH_ERR)
        return super().__add__(other=other)

    def __neg__(self) -> '_PolynomialNTTRepresentation':
        return _PolynomialNTTRepresentation(mod=self.mod, vals=tuple([-x for x in self._vals]))

    def __mul__(self, other: '_PolynomialNTTRepresentation') -> '_PolynomialNTTRepresentation':
        if other == 1:
            return self
        elif not _is_mul_valid(left=self, right=other):
            raise TypeError(_TYPE_MISMATCH_ERR)
        return _PolynomialNTTRepresentation(mod=self.mod, vals=tuple([cent(val=x * y, mod=self.mod) for x, y in zip(self._vals, other._vals)]))

    def __mod__(self, other: int) -> '_PolynomialNTTRepresentation':
        if not isinstance(other, int):
            raise TypeError(_MUST_BE_INT_ERR)
        return _PolynomialNTTRepresentation(mod=other, vals=tuple([cent(val=x, mod=other) for x in self._vals]))

    def transform(self) -> _PolynomialCoefficientRepresentation:
        return _PolynomialCoefficientRepresentation(mod=self.mod, vals=tuple(ntt(val=list(self._vals), mod=self.mod, inv_flag=True)))

    @property
    def coefs_norm_wght(self) -> Tuple[Tuple[int, ...], int, int]:
        return self.transform().coefs_norm_wght


class _PolynomialFactory:
    @staticmethod
    def make(mod: int, vals: Tuple[int, ...], rep_flag: str) -> Union[_PolynomialCoefficientRepresentation, _PolynomialNTTRepresentation]:
        if rep_flag == _POLYNOMIAL_REPRESENTATION_TYPES[0]:
            return _PolynomialCoefficientRepresentation(mod=mod, vals=vals)
        elif rep_flag == _POLYNOMIAL_REPRESENTATION_TYPES[1]:
            return _PolynomialNTTRepresentation(mod=mod, vals=vals)
        raise ValueError(_INVALID_REP_TYPE_ERR + f": {rep_flag}")


class Polynomial:
    ntt_rep: _PolynomialNTTRepresentation

    def __init__(self, mod: int, vals: Tuple[int, ...], rep_flag: str = _POLYNOMIAL_REPRESENTATION_TYPES[0]):
        """
        Initialize a polynomial with the given modulus, list of values (coefficients),
        and representation type ("coefficient" or "ntt").
        """
        if rep_flag not in _POLYNOMIAL_REPRESENTATION_TYPES:
            raise ValueError(_INVALID_REP_TYPE_ERR)
        elif rep_flag == _POLYNOMIAL_REPRESENTATION_TYPES[0]:
            x = _PolynomialFactory.make(mod=mod, vals=vals, rep_flag=rep_flag)
            self.ntt_rep = x.transform()
        elif rep_flag == _POLYNOMIAL_REPRESENTATION_TYPES[1]:
            self.ntt_rep = _PolynomialFactory.make(mod=mod, vals=vals, rep_flag=rep_flag)

    @property
    def mod(self) -> int:
        """
        Returns the modulus of the polynomial.
        """
        return self.ntt_rep.mod

    @property
    def vals(self) -> Tuple[int, ...]:
        """
        Returns the values (coefficients or NTT form) of the polynomial.
        """
        return self.ntt_rep._vals

    @property
    def deg(self) -> int:
        """
        Returns the degree of the polynomial.
        """
        return self.ntt_rep.deg

    @property
    def root_order(self) -> int:
        """
        Returns the root order of the polynomial.
        """
        return self.ntt_rep.root_order

    @property
    def root(self) -> int:
        """
        Returns the root of the polynomial.
        """
        return self.ntt_rep.root

    @property
    def inv_root(self) -> int:
        """
        Returns the inverse root of the polynomial.
        """
        return self.ntt_rep.inv_root

    @property
    def coef_norm_wght(self) -> Tuple[Tuple[int, ...], int, int]:
        coefs_rep: _PolynomialCoefficientRepresentation = self.ntt_rep.transform()
        return coefs_rep.coefs_norm_wght

    def __str__(self) -> str:
        """
        Returns a string representation of this polynomial.
        """
        return _POLYNOMIAL_REPRESENTATION_STR_PREFIX+f"(ntt={self.ntt_rep})"

    def __eq__(self, other) -> bool:
        """
        Checks if this polynomial is equal to another polynomial.
        """
        if not isinstance(other, Polynomial):
            return False
        return self.ntt_rep == other.ntt_rep

    def __add__(self, other: 'Polynomial') -> 'Polynomial':
        """
        Adds this polynomial with another polynomial, returning the result as a new Polynomial.
        """
        new_ntt_rep = self.ntt_rep + other.ntt_rep
        return Polynomial(mod=new_ntt_rep.mod, vals=new_ntt_rep._vals, rep_flag=_POLYNOMIAL_REPRESENTATION_TYPES[1])

    def __sub__(self, other: 'Polynomial') -> 'Polynomial':
        """
        Subtracts another polynomial from this polynomial, returning the result as a new Polynomial.
        """
        new_ntt_rep = self.ntt_rep - other.ntt_rep
        return Polynomial(mod=new_ntt_rep.mod, vals=new_ntt_rep._vals, rep_flag=_POLYNOMIAL_REPRESENTATION_TYPES[1])

    def __mul__(self, other: 'Polynomial') -> 'Polynomial':
        """
        Multiplies this polynomial with another polynomial, returning the result as a new Polynomial.
        """
        new_ntt_rep = self.ntt_rep * other.ntt_rep
        return Polynomial(mod=new_ntt_rep.mod, vals=new_ntt_rep._vals, rep_flag=_POLYNOMIAL_REPRESENTATION_TYPES[1])

    def __neg__(self) -> 'Polynomial':
        """
        Returns the negation of this polynomial.
        """
        new_ntt_rep = -self.ntt_rep
        return Polynomial(mod=new_ntt_rep.mod, vals=new_ntt_rep._vals, rep_flag=_POLYNOMIAL_REPRESENTATION_TYPES[1])

