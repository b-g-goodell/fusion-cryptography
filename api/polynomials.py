from typing import List as _List, Tuple as _Tuple
from algebra.polynomials import _PolynomialCoefficientRepresentation as _PolyC, _PolynomialNTTRepresentation as _PolyN, _PolynomialFactory as _Factory, _POLYNOMIAL_REPRESENTATION_TYPES
from api.errors import _INVALID_REP_TYPE_ERR


_POLYNOMIAL_REPRESENTATION_STR_PREFIX: str = 'Polynomial'


class Polynomial:
    ntt_rep: _PolyN

    def __init__(self, modulus: int, values: _List[int], representation: str = _POLYNOMIAL_REPRESENTATION_TYPES[0]):
        """
        Initialize a polynomial with the given modulus, list of values (coefficients),
        and representation type ("coefficient" or "ntt").
        """
        if representation not in _POLYNOMIAL_REPRESENTATION_TYPES:
            raise ValueError(_INVALID_REP_TYPE_ERR)
        elif representation == _POLYNOMIAL_REPRESENTATION_TYPES[0]:
            x = _Factory.create_representation(modulus=modulus, values=values, representation_type=representation)
            self.ntt_rep = x.transform()
        elif representation == _POLYNOMIAL_REPRESENTATION_TYPES[1]:
            self.ntt_rep = _Factory.create_representation(modulus=modulus, values=values, representation_type=representation)

    @property
    def modulus(self) -> int:
        """
        Returns the modulus of the polynomial.
        """
        return self.ntt_rep.modulus

    @property
    def values(self) -> _List[int]:
        """
        Returns the values (coefficients or NTT form) of the polynomial.
        """
        return self.ntt_rep.values

    @property
    def degree(self) -> int:
        """
        Returns the degree of the polynomial.
        """
        return self.ntt_rep.degree

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
    def coefs_norm_weight(self) -> _Tuple[_List[int], int, int]:
        coefs_rep: _PolyC = self.ntt_rep.transform()
        norm, weight = coefs_rep.norm_weight
        return coefs_rep.values, norm, weight

    def __str__(self):
        """
        Returns a string representation of this polynomial.
        """
        return _POLYNOMIAL_REPRESENTATION_STR_PREFIX+f"(ntt={self.ntt_rep})"

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other) -> bool:
        """
        Checks if this polynomial is equal to another polynomial.
        """
        if not isinstance(other, Polynomial):
            return False
        return self.ntt_rep == other.ntt_rep

    def __add__(self, other):
        """
        Adds this polynomial with another polynomial, returning the result as a new Polynomial.
        """
        new_ntt_rep = self.ntt_rep + other.ntt_rep
        return Polynomial(modulus=new_ntt_rep.modulus, values=new_ntt_rep.values, representation=_POLYNOMIAL_REPRESENTATION_TYPES[1])

    def __sub__(self, other):
        """
        Subtracts another polynomial from this polynomial, returning the result as a new Polynomial.
        """
        new_ntt_rep = self.ntt_rep - other.ntt_rep
        return Polynomial(modulus=new_ntt_rep.modulus, values=new_ntt_rep.values, representation=_POLYNOMIAL_REPRESENTATION_TYPES[1])

    def __mul__(self, other):
        """
        Multiplies this polynomial with another polynomial, returning the result as a new Polynomial.
        """
        new_ntt_rep = self.ntt_rep * other.ntt_rep
        return Polynomial(modulus=new_ntt_rep.modulus, values=new_ntt_rep.values, representation=_POLYNOMIAL_REPRESENTATION_TYPES[1])

    def __neg__(self):
        """
        Returns the negation of this polynomial.
        """
        new_ntt_rep = -self.ntt_rep
        return Polynomial(modulus=new_ntt_rep.modulus, values=new_ntt_rep.values, representation=_POLYNOMIAL_REPRESENTATION_TYPES[1])
