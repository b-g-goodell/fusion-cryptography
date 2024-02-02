# algebra/polynomials.py
from typing import List, Tuple, Union
from algebra.ntt import ntt, cent, is_ntt_friendly, find_prou
from algebra.errors import (_MUST_BE_INT_ERR, _MUST_HAVE_PROU_ERR, _MUST_BE_LIST_ERR, _TYPE_MISMATCH_ERR,
                            _MODULUS_MISMATCH_ERR, _NTT_NOT_IMPLEMENTED_ERR, _NORM_NOT_IMPLEMENTED_ERR,
                            _MUL_BASE_NOT_IMPLEMENTED_ERR, _INVALID_REP_TYPE_ERR, _DEGREE_MISMATCH_ERR, _MUST_BE_TUPLE_ERR)


# str prefixes
_POLYNOMIAL_REPRESENTATION_TYPES: List[str] = ['coefficient', 'ntt']
_POLYNOMIAL_REPRESENTATION_STR_PREFIX: str = '_PolynomialRepresentation'
_POLYNOMIAL_COEFFICIENT_REPRESENTATION_STR_PREFIX: str = '_PolynomialCoefficientRepresentation'
_POLYNOMIAL_NTT_REPRESENTATION_STR_PREFIX: str = '_PolynomialNTTRepresentation'
_POLYNOMIAL_REPRESENTATION_STR_PREFIX: str = 'Polynomial'


# types
TUPLE_OF_INTEGER_VALS: type = Tuple[int, ...]
TUPLE_OF_COEF_NORM_WGHT: type = Tuple[TUPLE_OF_INTEGER_VALS, int, int]


def _is_arithmetic_valid(left: '_PolynomialRepresentation', rght: '_PolynomialRepresentation') -> bool:
    return (type(left) is type(rght)
            and left.modulus == rght.modulus
            and left.deg == rght.deg
            and left.root_order == rght.root_order
            and left.root == rght.root
            and left.inv_root == rght.inv_root)


class _PolynomialRepresentation:
    _modulus: int
    _vals: TUPLE_OF_INTEGER_VALS

    def __init__(self, mod: int, vals: TUPLE_OF_INTEGER_VALS):
        if not isinstance(vals, tuple):
            raise TypeError(_MUST_BE_TUPLE_ERR)
        elif not isinstance(mod, int) or not all(isinstance(i, int) for i in vals):
            raise TypeError(_MUST_BE_INT_ERR)
        elif not is_ntt_friendly(mod=mod, deg=len(vals)):
            raise ValueError(_MUST_HAVE_PROU_ERR)
        self._modulus = mod
        self._vals = vals

    @property
    def modulus(self) -> int:
        return self._modulus

    @property
    def vals(self) -> TUPLE_OF_INTEGER_VALS:
        return self._vals

    @property
    def deg(self) -> int:
        return len(self._vals)

    @property
    def root_order(self) -> int:
        return 2*self.deg

    @property
    def root(self) -> int:
        return find_prou(mod=self.modulus, deg=self.deg)

    @property
    def inv_root(self) -> int:
        return pow(base=self.root, exp=self.modulus - 2, mod=self.modulus)

    @property
    def coef_norm_wght(self) -> TUPLE_OF_COEF_NORM_WGHT:
        raise NotImplementedError(_NORM_NOT_IMPLEMENTED_ERR)

    def __str__(self) -> str:
        return _POLYNOMIAL_REPRESENTATION_STR_PREFIX + f"(modulus={self.modulus}, values={self.vals})"

    def __eq__(self, other) -> bool:
        return (all((x - y) % self.modulus == 0 for x, y in zip(self.vals, other.vals))
                and _is_arithmetic_valid(left=self, rght=other))

    def __add__(self, other: '_PolynomialRepresentation') -> '_PolynomialRepresentation':
        if not _is_arithmetic_valid(left=self, rght=other):
            raise TypeError(_TYPE_MISMATCH_ERR)
        return _PolynomialRepresentation(mod=self.modulus, vals=tuple([cent(val=x + y, mod=self.modulus) for x, y in zip(self.vals, other.vals)]))

    def __radd__(self, other) -> '_PolynomialRepresentation':
        if other == 0:
            return self
        return self + other

    def __neg__(self) -> '_PolynomialRepresentation':
        return _PolynomialRepresentation(mod=self.modulus, vals=tuple([-x for x in self.vals]))

    def __sub__(self, other: '_PolynomialRepresentation') -> '_PolynomialRepresentation':
        return self + (-other)

    def __mul__(self, other: '_PolynomialRepresentation') -> '_PolynomialRepresentation':
        raise NotImplementedError(_MUL_BASE_NOT_IMPLEMENTED_ERR)

    def __mod__(self, other: int) -> '_PolynomialRepresentation':
        if not isinstance(other, int):
            raise TypeError(_MUST_BE_INT_ERR)
        return _PolynomialRepresentation(mod=other, vals=tuple([cent(val=x, mod=other) for x in self.vals]))


class _PolynomialCoefficientRepresentation(_PolynomialRepresentation):
    def __init__(self, mod: int, vals: TUPLE_OF_INTEGER_VALS):
        super().__init__(mod=mod, vals=vals)

    @property
    def coef_norm_wght(self) -> Tuple[TUPLE_OF_INTEGER_VALS, int, int]:
        return self.vals, max(abs(x) for x in self.vals), sum(1 if x % self.modulus != 0 else 0 for x in self.vals)

    def __str__(self) -> str:
        return _POLYNOMIAL_COEFFICIENT_REPRESENTATION_STR_PREFIX + f"(modulus={self.modulus}, degree={self.deg}, root={self.root}, inv_root={self.inv_root}, values={self.vals})"

    def __add__(self, other: '_PolynomialCoefficientRepresentation') -> '_PolynomialCoefficientRepresentation':
        if not isinstance(other, _PolynomialCoefficientRepresentation):
            raise TypeError(_TYPE_MISMATCH_ERR)
        x = super().__add__(other=other)
        return _PolynomialCoefficientRepresentation(mod=x.modulus, vals=x.vals)

    def __neg__(self) -> '_PolynomialCoefficientRepresentation':
        x = super().__neg__()
        return _PolynomialCoefficientRepresentation(mod=self.modulus, vals=x.vals)

    def __mul__(self, other: '_PolynomialCoefficientRepresentation') -> '_PolynomialCoefficientRepresentation':
        # TODO: multiply coef rep against ntt rep by transforming first
        if not isinstance(other, _PolynomialCoefficientRepresentation):
            raise TypeError(_TYPE_MISMATCH_ERR)
        elif not _is_arithmetic_valid(left=self, rght=other):
            raise TypeError(_TYPE_MISMATCH_ERR)
        new_values: List[int] = [0 for _ in range(2 * self.deg)]
        for i, x in enumerate(self.vals):
            for j, y in enumerate(other.vals):
                new_values[i + j] += x * y
        new_values = [cent(val=x - y, mod=self.modulus) for x, y in zip(new_values[:self.deg], new_values[self.deg:])]
        return _PolynomialCoefficientRepresentation(mod=self.modulus, vals=tuple(new_values))

    def __mod__(self, other: int) -> '_PolynomialCoefficientRepresentation':
        if not isinstance(other, int):
            raise TypeError(_MUST_BE_INT_ERR)
        return _PolynomialCoefficientRepresentation(mod=other, vals=tuple([cent(val=x, mod=other) for x in self.vals]))

    def transform_to_ntt_rep(self) -> '_PolynomialNTTRepresentation':
        return _PolynomialNTTRepresentation(mod=self.modulus, vals=tuple(ntt(val=list(self.vals), mod=self.modulus, inv_flag=False)))


class _PolynomialNTTRepresentation(_PolynomialRepresentation):
    def __init__(self, mod: int, vals: TUPLE_OF_INTEGER_VALS):
        super().__init__(mod=mod, vals=vals)

    @property
    def coef_norm_wght(self) -> Tuple[TUPLE_OF_INTEGER_VALS, int, int]:
        return self.transform_to_coef_rep().coef_norm_wght

    def __str__(self) -> str:
        return _POLYNOMIAL_NTT_REPRESENTATION_STR_PREFIX + f"(modulus={self.modulus}, degree={self.deg}, root={self.root}, inv_root={self.inv_root}, values={self.vals})"

    def __add__(self, other: '_PolynomialNTTRepresentation') -> '_PolynomialNTTRepresentation':
        if not isinstance(other, _PolynomialNTTRepresentation):
            raise TypeError(_TYPE_MISMATCH_ERR)
        x = super().__add__(other=other)
        return _PolynomialNTTRepresentation(mod=x.modulus, vals=x.vals)

    def __neg__(self) -> '_PolynomialNTTRepresentation':
        return _PolynomialNTTRepresentation(mod=self.modulus, vals=tuple([-x for x in self.vals]))

    def __mul__(self, other: '_PolynomialNTTRepresentation') -> '_PolynomialNTTRepresentation':
        if not isinstance(other, _PolynomialNTTRepresentation):
            raise TypeError(_TYPE_MISMATCH_ERR)
        elif not _is_arithmetic_valid(left=self, rght=other):
            raise TypeError(_TYPE_MISMATCH_ERR)
        return _PolynomialNTTRepresentation(mod=self.modulus, vals=tuple([cent(val=x * y, mod=self.modulus) for x, y in zip(self.vals, other.vals)]))

    def __mod__(self, other: int) -> '_PolynomialNTTRepresentation':
        if not isinstance(other, int):
            raise TypeError(_MUST_BE_INT_ERR)
        return _PolynomialNTTRepresentation(mod=other, vals=tuple([cent(val=x, mod=other) for x in self.vals]))

    def transform_to_coef_rep(self) -> _PolynomialCoefficientRepresentation:
        return _PolynomialCoefficientRepresentation(mod=self.modulus, vals=tuple(ntt(val=list(self.vals), mod=self.modulus, inv_flag=True)))


class _PolynomialFactory:
    @staticmethod
    def make(mod: int, vals: TUPLE_OF_INTEGER_VALS, rep_flag: str) -> Union[_PolynomialCoefficientRepresentation, _PolynomialNTTRepresentation]:
        if rep_flag == _POLYNOMIAL_REPRESENTATION_TYPES[0]:
            return _PolynomialCoefficientRepresentation(mod=mod, vals=vals)
        elif rep_flag == _POLYNOMIAL_REPRESENTATION_TYPES[1]:
            return _PolynomialNTTRepresentation(mod=mod, vals=vals)
        raise ValueError(_INVALID_REP_TYPE_ERR + f": {rep_flag}")


class Polynomial:
    ntt_rep: _PolynomialNTTRepresentation

    def __init__(self, mod: int, vals: TUPLE_OF_INTEGER_VALS, rep_flag: str = _POLYNOMIAL_REPRESENTATION_TYPES[0]):
        """
        Initialize a polynomial with the given modulus, list of values (coefficients),
        and representation type ("coefficient" or "ntt").
        """
        if rep_flag not in _POLYNOMIAL_REPRESENTATION_TYPES:
            raise ValueError(_INVALID_REP_TYPE_ERR)
        elif rep_flag == _POLYNOMIAL_REPRESENTATION_TYPES[0]:
            x = _PolynomialFactory.make(mod=mod, vals=vals, rep_flag=rep_flag)
            self.ntt_rep = x.transform_to_ntt_rep()
        elif rep_flag == _POLYNOMIAL_REPRESENTATION_TYPES[1]:
            self.ntt_rep = _PolynomialFactory.make(mod=mod, vals=vals, rep_flag=rep_flag)

    @property
    def modulus(self) -> int:
        """
        Returns the modulus of the polynomial.
        """
        return self.ntt_rep.modulus

    @property
    def vals(self) -> TUPLE_OF_INTEGER_VALS:
        """
        Returns the values (coefficients or NTT form) of the polynomial.
        """
        return self.ntt_rep._vals

    @property
    def deg(self) -> int:
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
    def coef_norm_wght(self) -> Tuple[TUPLE_OF_INTEGER_VALS, int, int]:
        coefs_rep: _PolynomialCoefficientRepresentation = self.ntt_rep.transform_to_coef_rep()
        return coefs_rep.coef_norm_wght

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
        return other.ntt_rep == self.ntt_rep

    def __add__(self, other: 'Polynomial') -> 'Polynomial':
        """
        Adds this polynomial with another polynomial, returning the result as a new Polynomial.
        """
        new_ntt_rep = self.ntt_rep + other.ntt_rep
        return Polynomial(mod=new_ntt_rep.modulus, vals=new_ntt_rep._vals, rep_flag=_POLYNOMIAL_REPRESENTATION_TYPES[1])

    def __sub__(self, other: 'Polynomial') -> 'Polynomial':
        """
        Subtracts another polynomial from this polynomial, returning the result as a new Polynomial.
        """
        new_ntt_rep = self.ntt_rep - other.ntt_rep
        return Polynomial(mod=new_ntt_rep.modulus, vals=new_ntt_rep._vals, rep_flag=_POLYNOMIAL_REPRESENTATION_TYPES[1])

    def __mul__(self, other: 'Polynomial') -> 'Polynomial':
        """
        Multiplies this polynomial with another polynomial, returning the result as a new Polynomial.
        """
        new_ntt_rep = self.ntt_rep * other.ntt_rep
        return Polynomial(mod=new_ntt_rep.modulus, vals=new_ntt_rep._vals, rep_flag=_POLYNOMIAL_REPRESENTATION_TYPES[1])

    def __neg__(self) -> 'Polynomial':
        """
        Returns the negation of this polynomial.
        """
        new_ntt_rep = -self.ntt_rep
        return Polynomial(mod=new_ntt_rep.modulus, vals=new_ntt_rep._vals, rep_flag=_POLYNOMIAL_REPRESENTATION_TYPES[1])
