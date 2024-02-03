"""
algebra/polynomials.py

A Polynomial class for handling fast arithmetic in constant-time.

Author: Brandon G Goodell
"""
from typing import List, Tuple, Union, Any, Type
from algebra.ntt import ntt, cent, is_ntt_friendly, find_prou
from algebra.errors import (_MUST_BE_INT_ERR, _MUST_HAVE_PROU_ERR, _TYPE_MISMATCH_ERR, _NORM_NOT_IMPLEMENTED_ERR,
                            _MUL_BASE_NOT_IMPLEMENTED_ERR, _INVALID_REP_TYPE_ERR, _MUST_BE_TUPLE_ERR)
import struct


# allowable representations
_ALLOWED_REPS: List[str] = ['coefficient', 'ntt']

# byte order for serialization
BIG_ENDIAN: str = '>'


# __str__ prefixes
_POLYNOMIAL_REPRESENTATION_STR_PREFIX: str = '_PolynomialRepresentation'
_POLYNOMIAL_COEFFICIENT_REPRESENTATION_STR_PREFIX: str = '_PolynomialCoefficientRepresentation'
_POLYNOMIAL_NTT_REPRESENTATION_STR_PREFIX: str = '_PolynomialNTTRepresentation'


# types
IntegerTuple: type = Tuple[int, ...]
CoefNormWghtTuple: type = Tuple[IntegerTuple, int, int]


def _is_same_type(left: Any, rght: Any) -> bool:
    """
    Returns true if left and rght are of the same type, false otherwise.

    :param left: The left-hand side of the comparison.
    :type left: Any
    :param rght: The right-hand side of the comparison.
    :type rght: Any
    :returns: True if left and rght are of the same type, false otherwise.
    :rtype: bool
    """
    return type(left) is type(rght)


def _is_same_ring(left: '_PolynomialRepresentation', rght: '_PolynomialRepresentation') -> bool:
    """
    Returns true if left and rght represent polynomials in the same ring, false otherwise.

    :param left: The left-hand side of the comparison.
    :type left: _PolynomialRepresentation
    :param rght: The right-hand side of the comparison.
    :type rght: _PolynomialRepresentation
    :returns: True if left and rght represent polynomials in the same ring, false otherwise.
    :rtype: bool
    """
    return (left.modulus == rght.modulus
            and left.deg == rght.deg
            and left.deg == len(left.vals)
            and rght.deg == len(rght.vals)
            and left.root_order == rght.root_order
            and left.root == rght.root
            and left.inv_root == rght.inv_root)


class _PolynomialRepresentation:
    """
    Base class for polynomial representations. Intended to be immutable and constant-time. Implements basic
    arithmetic operations and string representation.

    :ivar _modulus: The modulus of the polynomial.
    :type _modulus: int
    :ivar _vals: Representation of the polynomial. Interpretation of this representation varies (see child classes).
    :type _vals: IntegerTuple
    """
    
    _modulus: int
    _vals: IntegerTuple

    def __init__(self, mod: int, vals: IntegerTuple):
        if not isinstance(mod, int):
            raise TypeError(_MUST_BE_INT_ERR)
        elif not isinstance(vals, tuple):
            raise TypeError(_MUST_BE_TUPLE_ERR)
        elif not all(isinstance(i, int) for i in vals):
            raise TypeError(_MUST_BE_INT_ERR)
        elif not is_ntt_friendly(mod=mod, deg=len(vals)):
            raise ValueError(_MUST_HAVE_PROU_ERR)
        self._modulus = mod
        self._vals = vals

    @property
    def modulus(self) -> int:
        ''' Returns the modulus in the ideal (modulus, X**d + 1).'''
        return self._modulus

    @property
    def vals(self) -> IntegerTuple:
        '''
        Returns the numerical representation of the polynomial. Interpretation of this representation varies
        (see child classes).
        '''
        return self._vals

    @property
    def deg(self) -> int:
        ''' Returns the degree d of the polynomial generator of the ideal (modulus, X**d + 1).'''
        return len(self._vals)

    @property
    def root_order(self) -> int:
        ''' Returns the order of the primitive root of unity used for NTT and INTT transforms.'''
        return 2*self.deg

    @property
    def root(self) -> int:
        ''' Returns the primitive root of unity used for NTT and INTT transforms. '''
        return find_prou(mod=self.modulus, deg=self.deg)

    @property
    def inv_root(self) -> int:
        ''' Returns the inverse of the primitive root of unity used for NTT and INTT transforms. '''
        return pow(base=self.root, exp=self.modulus - 2, mod=self.modulus)

    @property
    def coef_norm_wght(self) -> CoefNormWghtTuple:
        ''' Returns the coefficients, norm, and weight of the polynomial.'''
        raise NotImplementedError(_NORM_NOT_IMPLEMENTED_ERR)

    def __str__(self) -> str:
        return _POLYNOMIAL_REPRESENTATION_STR_PREFIX + f"(modulus={self.modulus}, values={self.vals})"

    def __eq__(self, other) -> bool:
        return (_is_same_type(left=self, rght=other) and _is_same_ring(left=self, rght=other)
                and all((x - y) % self.modulus == 0 for x, y in zip(self.vals, other.vals)))

    def __add__(self, other: '_PolynomialRepresentation') -> '_PolynomialRepresentation':
        if not _is_same_type(left=self, rght=other):
            raise TypeError(_TYPE_MISMATCH_ERR)
        elif not _is_same_ring(left=self, rght=other):
            raise TypeError(_TYPE_MISMATCH_ERR)
        return type(self)(mod=self.modulus, vals=tuple(cent(val=x + y, mod=self.modulus) for x, y in zip(self.vals, other.vals)))

    def __radd__(self, other) -> '_PolynomialRepresentation':
        ''' Intended to allow usage of sum() but does not appear to be working? Defunct this? '''
        if other == 0:
            return self
        return self + other

    def __neg__(self) -> '_PolynomialRepresentation':

        return _PolynomialRepresentation(mod=self.modulus, vals=tuple([-x for x in self.vals]))

    def __sub__(self, other: '_PolynomialRepresentation') -> '_PolynomialRepresentation':
        return self + (-other)

    def __mul__(self, other: '_PolynomialRepresentation') -> '_PolynomialRepresentation':
        if not _is_same_type(left=self, rght=other):
            raise TypeError(_TYPE_MISMATCH_ERR)
        elif not _is_same_ring(left=self, rght=other):
            raise TypeError(_TYPE_MISMATCH_ERR)
        raise NotImplementedError(_MUL_BASE_NOT_IMPLEMENTED_ERR)

    def __mod__(self, other: int) -> '_PolynomialRepresentation':
        '''
        Mod out a polynomial by a new modulus. Returns self when the new modulus is greater than or equal to the
        current modulus.
        '''
        if not isinstance(other, int):
            raise TypeError(_MUST_BE_INT_ERR)
        elif other > self.modulus//2:
            raise ValueError("New modulus must be less than half the current modulus for __mod__ function to be interpreted as surjection.")
        result_vals: Tuple[int, ...] = tuple(cent(val=x, mod=other) for x in self.vals)
        return _PolynomialRepresentation(mod=other, vals=result_vals)


    def to_bytes(self) -> bytes:
        '''
        Serializes the polynomial representation to bytes.
        :returns: The bytes representation of the polynomial.
        :rtype: bytes
        '''
        d: int = self.deg.bit_length() - 1
        mod_bytes: bytes = struct.pack(BIG_ENDIAN + 'I', self.modulus)
        d_byte: bytes = struct.pack('B', d)
        vals_bytes: bytes = b''.join(struct.pack(BIG_ENDIAN + 'i', val) for val in self.vals)
        return mod_bytes + d_byte + vals_bytes

    @classmethod
    def from_bytes(cls, data: bytes) -> '_PolynomialRepresentation':
        '''
        Deserializes a polynomial representation from bytes.
        :param data: The bytes representation of the polynomial.
        :type data: bytes
        :returns: A polynomial representation instance.
        :rtype: _PolynomialRepresentation
        :raises ValueError: If the bytes representation is invalid.
        '''
        if len(data) < 5 or (len(data) - 5) % 4 != 0:
            raise ValueError("Invalid bytes representation")
        mod: int = struct.unpack(BIG_ENDIAN + 'I', data[:4])[0]
        d: int = struct.unpack('B', data[4:5])[0]
        len_vals: int = 1 << d
        if len(data) - 5 != 4 * len_vals:
            raise ValueError("Invalid bytes representation")
        vals: Tuple[int, ...] = struct.unpack(BIG_ENDIAN + 'i' * len_vals, data[5:])  # Remaining bytes for values
        return cls(mod=mod, vals=vals)


class _PolynomialCoefficientRepresentation(_PolynomialRepresentation):
    '''
    Represent a polynomial using a tuple of coefficients. Intended to be immutable and constant-time. Implements
    "chalkboard" or "FOIL" arithmetic for polynomial multiplication.
    '''
    def __init__(self, mod: int, vals: IntegerTuple):
        super().__init__(mod=mod, vals=vals)

    @property
    def coef_norm_wght(self) -> Tuple[IntegerTuple, int, int]:
        return self.vals, max(abs(x) for x in self.vals), sum(1 if x % self.modulus != 0 else 0 for x in self.vals)

    def __str__(self) -> str:
        return _POLYNOMIAL_COEFFICIENT_REPRESENTATION_STR_PREFIX + f"(modulus={self.modulus}, degree={self.deg}, root={self.root}, inv_root={self.inv_root}, values={self.vals})"

    def __add__(self, other: '_PolynomialCoefficientRepresentation') -> '_PolynomialCoefficientRepresentation':
        if not _is_same_type(left=self, rght=other):
            raise TypeError(_TYPE_MISMATCH_ERR)
        elif not _is_same_ring(left=self, rght=other):
            raise TypeError(_TYPE_MISMATCH_ERR)
        x = super().__add__(other=other)
        return _PolynomialCoefficientRepresentation(mod=x.modulus, vals=x.vals)

    def __neg__(self) -> '_PolynomialCoefficientRepresentation':
        x = super().__neg__()
        return _PolynomialCoefficientRepresentation(mod=self.modulus, vals=x.vals)

    def __mul__(self, other: '_PolynomialCoefficientRepresentation') -> '_PolynomialCoefficientRepresentation':
        if not _is_same_type(left=self, rght=other):
            raise TypeError(_TYPE_MISMATCH_ERR)
        elif not _is_same_ring(left=self, rght=other):
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
        elif other >= self.modulus//2:
            return self
        return _PolynomialCoefficientRepresentation(mod=other, vals=tuple([cent(val=x, mod=other) for x in self.vals]))

    def transform_to_ntt_rep(self) -> '_PolynomialNTTRepresentation':
        return _PolynomialNTTRepresentation(mod=self.modulus, vals=tuple(ntt(val=list(self.vals), mod=self.modulus, inv_flag=False)))


class _PolynomialNTTRepresentation(_PolynomialRepresentation):
    def __init__(self, mod: int, vals: IntegerTuple):
        super().__init__(mod=mod, vals=vals)

    @property
    def coef_norm_wght(self) -> Tuple[IntegerTuple, int, int]:
        return self.transform_to_coef_rep().coef_norm_wght

    def __str__(self) -> str:
        return _POLYNOMIAL_NTT_REPRESENTATION_STR_PREFIX + f"(modulus={self.modulus}, degree={self.deg}, root={self.root}, inv_root={self.inv_root}, values={self.vals})"

    def __add__(self, other: '_PolynomialNTTRepresentation') -> '_PolynomialNTTRepresentation':
        if not _is_same_type(left=self, rght=other):
            raise TypeError(_TYPE_MISMATCH_ERR)
        elif not _is_same_ring(left=self, rght=other):
            raise TypeError(_TYPE_MISMATCH_ERR)
        x = super().__add__(other=other)
        return _PolynomialNTTRepresentation(mod=x.modulus, vals=x.vals)

    def __neg__(self) -> '_PolynomialNTTRepresentation':
        return _PolynomialNTTRepresentation(mod=self.modulus, vals=tuple([-x for x in self.vals]))

    def __mul__(self, other: '_PolynomialNTTRepresentation') -> '_PolynomialNTTRepresentation':
        if not _is_same_type(left=self, rght=other):
            raise TypeError(_TYPE_MISMATCH_ERR)
        elif not _is_same_ring(left=self, rght=other):
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
    def make(mod: int, vals: IntegerTuple, rep_flag: str) -> Union[_PolynomialCoefficientRepresentation, _PolynomialNTTRepresentation]:
        if rep_flag == _ALLOWED_REPS[0]:
            return _PolynomialCoefficientRepresentation(mod=mod, vals=vals)
        elif rep_flag == _ALLOWED_REPS[1]:
            return _PolynomialNTTRepresentation(mod=mod, vals=vals)
        raise ValueError(_INVALID_REP_TYPE_ERR + f": {rep_flag}")


class Polynomial:
    ntt_rep: _PolynomialNTTRepresentation

    def __init__(self, mod: int, vals: IntegerTuple, rep_flag: str = _ALLOWED_REPS[0]):
        """
        Initialize a polynomial with the given modulus, list of values (coefficients),
        and representation type ("coefficient" or "ntt").
        """
        if rep_flag not in _ALLOWED_REPS:
            raise ValueError(_INVALID_REP_TYPE_ERR)
        elif rep_flag == _ALLOWED_REPS[0]:
            x = _PolynomialFactory.make(mod=mod, vals=vals, rep_flag=rep_flag)
            self.ntt_rep = x.transform_to_ntt_rep()
        elif rep_flag == _ALLOWED_REPS[1]:
            self.ntt_rep = _PolynomialFactory.make(mod=mod, vals=vals, rep_flag=rep_flag)

    @property
    def modulus(self) -> int:
        """
        Returns the modulus of the polynomial.
        """
        return self.ntt_rep.modulus

    @property
    def vals(self) -> IntegerTuple:
        """
        Returns the values (coefficients or NTT form) of the polynomial.
        """
        return self.ntt_rep.vals

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
    def coef_norm_wght(self) -> Tuple[IntegerTuple, int, int]:
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
        return Polynomial(mod=new_ntt_rep.modulus, vals=new_ntt_rep.vals, rep_flag=_ALLOWED_REPS[1])

    def __sub__(self, other: 'Polynomial') -> 'Polynomial':
        """
        Subtracts another polynomial from this polynomial, returning the result as a new Polynomial.
        """
        new_ntt_rep = self.ntt_rep - other.ntt_rep
        return Polynomial(mod=new_ntt_rep.modulus, vals=new_ntt_rep.vals, rep_flag=_ALLOWED_REPS[1])

    def __mul__(self, other: 'Polynomial') -> 'Polynomial':
        """
        Multiplies this polynomial with another polynomial, returning the result as a new Polynomial.
        """
        new_ntt_rep = self.ntt_rep * other.ntt_rep
        return Polynomial(mod=new_ntt_rep.modulus, vals=new_ntt_rep.vals, rep_flag=_ALLOWED_REPS[1])

    def __neg__(self) -> 'Polynomial':
        """
        Returns the negation of this polynomial.
        """
        new_ntt_rep = -self.ntt_rep
        return Polynomial(mod=new_ntt_rep.modulus, vals=new_ntt_rep._vals, rep_flag=_ALLOWED_REPS[1])




