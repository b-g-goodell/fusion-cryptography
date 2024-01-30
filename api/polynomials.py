# api/polynomials.py
from typing import Tuple, Union
from algebra.polynomials import (PolynomialCoefficientRepresentation as PolyC, PolynomialNTTRepresentation as PolyN,
                                 _PolynomialFactory as _Factory, _POLYNOMIAL_REPRESENTATION_TYPES)
from api.errors import _INVALID_REP_TYPE_ERR


_POLYNOMIAL_REPRESENTATION_STR_PREFIX: str = 'Polynomial'


class Polynomial:
    ntt_rep: PolyN

    def __init__(self, mod: int, vals: Tuple[int, ...], rep_flag: str = _POLYNOMIAL_REPRESENTATION_TYPES[0]):
        """
        Initialize a polynomial with the given modulus, list of values (coefficients),
        and representation type ("coefficient" or "ntt").
        """
        if rep_flag not in _POLYNOMIAL_REPRESENTATION_TYPES:
            raise ValueError(_INVALID_REP_TYPE_ERR)
        elif rep_flag == _POLYNOMIAL_REPRESENTATION_TYPES[0]:
            x = _Factory.make(mod=mod, vals=vals, rep_flag=rep_flag)
            self.ntt_rep = x.transform()
        elif rep_flag == _POLYNOMIAL_REPRESENTATION_TYPES[1]:
            self.ntt_rep = _Factory.make(mod=mod, vals=vals, rep_flag=rep_flag)

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
        return self.ntt_rep.vals

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
    def coefs_norm_weight(self) -> Tuple[Tuple[int, ...], int, int]:
        coefs_rep: PolyC = self.ntt_rep.transform()
        return coefs_rep.coefs_norm_wght

    def __str__(self) -> str:
        """
        Returns a string representation of this polynomial.
        """
        return _POLYNOMIAL_REPRESENTATION_STR_PREFIX+f"(ntt={self.ntt_rep})"

    def __repr__(self) -> str:
        return self.__str__()

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
        return Polynomial(mod=new_ntt_rep.mod, vals=new_ntt_rep.vals, rep_flag=_POLYNOMIAL_REPRESENTATION_TYPES[1])

    def __sub__(self, other: 'Polynomial') -> 'Polynomial':
        """
        Subtracts another polynomial from this polynomial, returning the result as a new Polynomial.
        """
        new_ntt_rep = self.ntt_rep - other.ntt_rep
        return Polynomial(mod=new_ntt_rep.mod, vals=new_ntt_rep.vals, rep_flag=_POLYNOMIAL_REPRESENTATION_TYPES[1])

    def __mul__(self, other: 'Polynomial') -> 'Polynomial':
        """
        Multiplies this polynomial with another polynomial, returning the result as a new Polynomial.
        """
        new_ntt_rep = self.ntt_rep * other.ntt_rep
        return Polynomial(mod=new_ntt_rep.mod, vals=new_ntt_rep.vals, rep_flag=_POLYNOMIAL_REPRESENTATION_TYPES[1])

    def __neg__(self) -> 'Polynomial':
        """
        Returns the negation of this polynomial.
        """
        new_ntt_rep = -self.ntt_rep
        return Polynomial(mod=new_ntt_rep.mod, vals=new_ntt_rep.vals, rep_flag=_POLYNOMIAL_REPRESENTATION_TYPES[1])

