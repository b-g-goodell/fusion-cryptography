"""
algebra/polynomials.py

Handle polynomial arithmetic in constant time with polynomials from Z_p[X]/(X^d + 1) for p a prime and d
a positive power of 2 such that (p-1) % (2d) == 0.

Author: Brandon G Goodell
"""
# TODO: add ser/deser
from algebra.errors import (_TypeMismatchError, _MustHavePROUError, _NormNotImplementedError, _MulBaseNotImplementedError, _MustBeIntError, _InvalidRepTypeError, _MustBeTupleError, _DecodingError)
from algebra.ntt import _is_ntt_friendly, _find_prou, _cent, _ntt

ENDIANNESS: str = 'big'

_POLYNOMIAL_REP_TYPES: list[str] = ['coefficient', 'ntt']


def _is_arithmetic_valid(left: '_PolynomialRep', rght: '_PolynomialRep') -> bool:
    """
    Checks if the arithmetic operations are valid between two polynomials.
    When computing f*g, we think of f as the left polynomial and g as the right.
    :param left: The left polynomial
    :type left: _PolynomialRep
    :param rght: The right polynomial
    :type rght: _PolynomialRep
    :return: Boolean indicating whether arithmetic is valid.
    :rtype: bool
    """
    return (type(left) is type(rght)
            and left.modulus == rght.modulus
            and left.deg == rght.deg)


class _PolynomialRep:
    """
    The base class for polynomial representations. This class should not be instantiated directly. It implements the
    modulus, vals, deg, root_order, root, inv_root, and coef_norm_wght properties, as well as the __str__, __repr__,
    __eq__, __add__, __radd__, __neg__, __sub__, and __mod__ methods. The __mul__ implementations
    vary by subclass.

    :ivar _modulus: The modulus of the polynomial.
    :type _modulus: int
    :ivar _vals: The values (coefficients or NTT form) of the polynomial.
    :type _vals: tuple[int, ...]
    """
    _modulus: int
    _vals: tuple[int, ...]

    def __init__(self, mod: int, vals: tuple[int, ...]):
        if not isinstance(vals, tuple):
            raise _MustBeTupleError
        elif not all(isinstance(i, int) for i in vals) or not isinstance(mod, int):
            raise _MustBeIntError
        elif not _is_ntt_friendly(mod=mod, deg=len(vals)):
            raise _MustHavePROUError
        self._modulus = mod
        self._vals = vals

    @property
    def modulus(self) -> int:
        return self._modulus

    @property
    def vals(self) -> tuple[int, ...]:
        return self._vals

    @property
    def deg(self) -> int:
        return len(self._vals)

    @property
    def root_order(self) -> int:
        return 2*self.deg

    @property
    def root(self) -> int:
        return _find_prou(mod=self.modulus, deg=self.deg)

    @property
    def inv_root(self) -> int:
        return pow(base=self.root, exp=self.modulus - 2, mod=self.modulus)

    @property
    def coef_norm_wght(self) -> tuple[tuple[int, ...], int, int]:
        raise _NormNotImplementedError

    def __str__(self) -> str:
        return self.__class__.__name__ + f"(modulus={self.modulus}, values={self.vals})"

    def __repr__(self) -> str:
        return self.__str__()

    def __eq__(self, other) -> bool:
        return (_is_arithmetic_valid(left=self, rght=other) and
                all((x - y) % self.modulus == 0 for x, y in zip(self.vals, other.vals)))

    def __add__(self, other: '_PolynomialRep') -> '_PolynomialRep':
        if not _is_arithmetic_valid(left=self, rght=other):
            raise _TypeMismatchError
        return _PolynomialRep(
            mod=self.modulus, vals=tuple([_cent(val=x + y, mod=self.modulus) for x, y in zip(self.vals, other.vals)]))

    def __radd__(self, other) -> '_PolynomialRep':
        if other == 0:
            return self
        return self + other

    def __neg__(self) -> '_PolynomialRep':
        return _PolynomialRep(mod=self.modulus, vals=tuple([-x for x in self.vals]))

    def __sub__(self, other: '_PolynomialRep') -> '_PolynomialRep':
        return self + (-other)

    def __mul__(self, other: '_PolynomialRep') -> '_PolynomialRep':
        raise _MulBaseNotImplementedError

    def __mod__(self, other: int) -> '_PolynomialRep':
        if not isinstance(other, int):
            raise _MustBeIntError
        return _PolynomialRep(mod=other, vals=tuple([_cent(val=x, mod=other) for x in self.vals]))

    def to_bytes(self) -> bytes:
        x: bytes = self.modulus.to_bytes(length=4, byteorder=ENDIANNESS) + len(self.vals).to_bytes(length=4, byteorder=ENDIANNESS)
        x += b''.join(i.to_bytes(length=4, byteorder=ENDIANNESS) for i in self.vals)
        c, n, w = self.coef_norm_wght
        x += n.to_bytes(length=4, byteorder=ENDIANNESS) + w.to_bytes(length=4, byteorder=ENDIANNESS)
        return x

    @staticmethod
    def from_bytes(b: bytes) -> '_PolynomialRep':
        mod: int = int.from_bytes(b[:4], byteorder=ENDIANNESS)
        deg: int = int.from_bytes(b[4:8], byteorder=ENDIANNESS)
        vals: tuple[int, ...] = tuple(int.from_bytes(b[i:i+4], byteorder=ENDIANNESS) for i in range(8, 8 + 4*deg))
        decoded_n: int = int.from_bytes(b[8 + 4*deg:12 + 4*deg], byteorder=ENDIANNESS)
        decoded_w: int = int.from_bytes(b[12 + 4*deg:16 + 4*deg], byteorder=ENDIANNESS)
        result: _PolynomialRep = _PolynomialRep(mod=mod, vals=vals)
        _, n, w = result.coef_norm_wght
        if n != decoded_n or w != decoded_w:
            raise _DecodingError
        return _PolynomialRep(mod=mod, vals=vals)


class _PolynomialCoefRep(_PolynomialRep):
    """
    The coefficient representation of a polynomial. This class can be instantiated directly, but multiplication
    is inefficient. This class wraps the __add__ and __neg__ classes from _PolynomialRepresentation
    with a simple type check, and implements multiplication by "FOIL"ing the coefficients.

    Also implements the transform_to_ntt_rep method, which returns a new _PolynomialNTTRep object.
    """
    def __init__(self, mod: int, vals: tuple[int, ...]):
        super().__init__(mod=mod, vals=vals)

    @property
    def coef_norm_wght(self) -> tuple[tuple[int, ...], int, int]:
        return self.vals, max(abs(x) for x in self.vals), sum(1 if x % self.modulus != 0 else 0 for x in self.vals)

    def __add__(self, other: '_PolynomialCoefRep') -> '_PolynomialCoefRep':
        if not isinstance(other, _PolynomialCoefRep):
            raise _TypeMismatchError
        x = super().__add__(other=other)
        return _PolynomialCoefRep(mod=x.modulus, vals=x.vals)

    def __neg__(self) -> '_PolynomialCoefRep':
        x = super().__neg__()
        return _PolynomialCoefRep(mod=self.modulus, vals=x.vals)

    def __mul__(self, other: '_PolynomialCoefRep') -> '_PolynomialCoefRep':
        # TODO: multiply coef rep against ntt rep by transforming first
        if not isinstance(other, _PolynomialCoefRep):
            raise _TypeMismatchError
        elif not _is_arithmetic_valid(left=self, rght=other):
            raise _TypeMismatchError
        new_values: list[int] = [0 for _ in range(2 * self.deg)]
        for i, x in enumerate(self.vals):
            for j, y in enumerate(other.vals):
                new_values[i + j] += x * y
        new_values = [_cent(val=x - y, mod=self.modulus) for x, y in zip(new_values[:self.deg], new_values[self.deg:])]
        return _PolynomialCoefRep(mod=self.modulus, vals=tuple(new_values))

    def transform_to_ntt_rep(self) -> '_PolynomialNTTRep':
        return _PolynomialNTTRep(
            mod=self.modulus, vals=tuple(_ntt(val=list(self.vals), mod=self.modulus, inv_flag=False)))


class _PolynomialNTTRep(_PolynomialRep):
    """
    The NTT representation of a polynomial. This class can be instantiated directly, and is more efficient for
    multiplication. This class wraps the __add__ and __neg__ classes from _PolynomialRepresentation
    with a simple type check, and implements multiplication by element-wise multiplication of the coefficients.
    Computing the coefficients, norm, and weight requires transforming back to the coefficient representation.

    This class also implements the transform_to_coef_rep method, which returns a new _PolynomialCoefRep object.
    """
    def __init__(self, mod: int, vals: tuple[int, ...]):
        super().__init__(mod=mod, vals=vals)

    @property
    def coef_norm_wght(self) -> tuple[tuple[int, ...], int, int]:
        return self.transform_to_coef_rep().coef_norm_wght

    def __add__(self, other: '_PolynomialNTTRep') -> '_PolynomialNTTRep':
        if not isinstance(other, _PolynomialNTTRep):
            raise _TypeMismatchError
        x = super().__add__(other=other)
        return _PolynomialNTTRep(mod=x.modulus, vals=x.vals)

    def __neg__(self) -> '_PolynomialNTTRep':
        return _PolynomialNTTRep(mod=self.modulus, vals=tuple([-x for x in self.vals]))

    def __mul__(self, other: '_PolynomialNTTRep') -> '_PolynomialNTTRep':
        if not isinstance(other, _PolynomialNTTRep):
            raise _TypeMismatchError
        elif not _is_arithmetic_valid(left=self, rght=other):
            raise _TypeMismatchError
        return _PolynomialNTTRep(
            mod=self.modulus, vals=tuple([_cent(val=x * y, mod=self.modulus) for x, y in zip(self.vals, other.vals)]))

    def __rmul__(self, other):
        if other == 1:
            return self
        return self * other

    def __mod__(self, other: int) -> '_PolynomialNTTRep':
        if not isinstance(other, int):
            raise _TypeMismatchError
        return _PolynomialNTTRep(mod=other, vals=tuple([_cent(val=x, mod=other) for x in self.vals]))

    def transform_to_coef_rep(self) -> _PolynomialCoefRep:
        return _PolynomialCoefRep(
            mod=self.modulus, vals=tuple(_ntt(val=list(self.vals), mod=self.modulus, inv_flag=True)))


class _PolynomialFactory:
    """
    A factory class for creating polynomials. This class should not be instantiated directly. It implements the _make
    method, which returns a new _PolynomialCoefRep or _PolynomialNTTRep object, depending on the representation type.
    """
    @staticmethod
    def _make(mod: int, vals: tuple[int, ...], rep_flag: str) -> _PolynomialCoefRep | _PolynomialNTTRep:
        if rep_flag == _POLYNOMIAL_REP_TYPES[0]:
            return _PolynomialCoefRep(mod=mod, vals=vals)
        elif rep_flag == _POLYNOMIAL_REP_TYPES[1]:
            return _PolynomialNTTRep(mod=mod, vals=vals)
        raise _InvalidRepTypeError


class _Polynomial:
    """
    A class for polynomials. This class should be instantiated directly. It can be instantiated with a coefficient
    representation or an NTT representation of a polynomial, automatically converting coefficient
    representations to NTT representations before storage. Under the hood, arithmetic efficiently takes place with the NTT
    representations while checking the coefficients, norm, and weights uses the coefficient representation.
    """
    _ntt_rep: _PolynomialNTTRep

    def __init__(self, mod: int, vals: tuple[int, ...], rep_flag: str = _POLYNOMIAL_REP_TYPES[0]):
        """
        Initialize a polynomial with the given modulus, list of values (coefficients),
        and representation type ("coefficient" or "ntt").
        """
        if rep_flag not in _POLYNOMIAL_REP_TYPES:
            raise _InvalidRepTypeError
        elif rep_flag == _POLYNOMIAL_REP_TYPES[0]:
            x = _PolynomialFactory._make(mod=mod, vals=vals, rep_flag=rep_flag)
            self._ntt_rep = x.transform_to_ntt_rep()
        elif rep_flag == _POLYNOMIAL_REP_TYPES[1]:
            self._ntt_rep = _PolynomialFactory._make(mod=mod, vals=vals, rep_flag=rep_flag)

    @property
    def modulus(self) -> int:
        """
        Returns the modulus of the polynomial.
        """
        return self._ntt_rep.modulus

    @property
    def vals(self) -> tuple[int, ...]:
        """
        Returns the values (coefficients or NTT form) of the polynomial.
        """
        return self._ntt_rep.vals

    @property
    def deg(self) -> int:
        return self._ntt_rep.deg

    @property
    def root_order(self) -> int:
        """
        Returns the root order of the polynomial.
        """
        return self._ntt_rep.root_order

    @property
    def root(self) -> int:
        """
        Returns the root of the polynomial.
        """
        return self._ntt_rep.root

    @property
    def inv_root(self) -> int:
        """
        Returns the inverse root of the polynomial.
        """
        return self._ntt_rep.inv_root

    @property
    def coef_norm_wght(self) -> tuple[tuple[int, ...], int, int]:
        coefs_rep: _PolynomialCoefRep = self._ntt_rep.transform_to_coef_rep()
        return coefs_rep.coef_norm_wght

    def __eq__(self, other) -> bool:
        """
        Checks if this polynomial is equal to another polynomial.
        """
        if not isinstance(other, _Polynomial):
            return False
        return other._ntt_rep == self._ntt_rep

    def __str__(self) -> str:
        return self.__class__.__name__ + f"(ntt={self._ntt_rep})"

    def __add__(self, other: '_Polynomial') -> '_Polynomial':
        """
        Adds this polynomial with another polynomial, returning the result as a new Polynomial.
        """
        new_ntt_rep = self._ntt_rep + other._ntt_rep
        return _Polynomial(mod=new_ntt_rep.modulus, vals=new_ntt_rep.vals, rep_flag=_POLYNOMIAL_REP_TYPES[1])

    def __sub__(self, other: '_Polynomial') -> '_Polynomial':
        """
        Subtracts another polynomial from this polynomial, returning the result as a new Polynomial.
        """
        new_ntt_rep = self._ntt_rep - other._ntt_rep
        return _Polynomial(mod=new_ntt_rep.modulus, vals=new_ntt_rep.vals, rep_flag=_POLYNOMIAL_REP_TYPES[1])

    def __mul__(self, other: '_Polynomial') -> '_Polynomial':
        """
        Multiplies this polynomial with another polynomial, returning the result as a new Polynomial.
        """
        new_ntt_rep = self._ntt_rep * other._ntt_rep
        return _Polynomial(mod=new_ntt_rep.modulus, vals=new_ntt_rep.vals, rep_flag=_POLYNOMIAL_REP_TYPES[1])

    def __neg__(self) -> '_Polynomial':
        """
        Returns the negation of this polynomial.
        """
        new_ntt_rep = -self._ntt_rep
        return _Polynomial(mod=new_ntt_rep.modulus, vals=new_ntt_rep.vals, rep_flag=_POLYNOMIAL_REP_TYPES[1])

    def to_bytes(self) -> bytes:
        return self._ntt_rep.to_bytes()

    @staticmethod
    def from_bytes_wo_norm_and_wght(b: bytes) -> '_Polynomial':
        if len(b) < 8:
            raise _DecodingError
        mod_bytes: bytes = b[:4]
        mod: int = int.from_bytes(mod_bytes, byteorder=ENDIANNESS)
        deg_bytes: bytes = b[4:8]
        deg: int = int.from_bytes(deg_bytes, byteorder=ENDIANNESS)
        if len(b) < 8 + 4*deg:
            raise _DecodingError
        vals_bytes: bytes = b[8:8 + 4*deg]
        vals: tuple[int, ...] = tuple(int.from_bytes(vals_bytes[i:i+4], byteorder=ENDIANNESS) for i in range(8, 8 + 4*deg))
        return _Polynomial(mod=mod, vals=vals, rep_flag=_POLYNOMIAL_REP_TYPES[1])

    @staticmethod
    def from_bytes(b: bytes) -> '_Polynomial':
        decoded_poly: _Polynomial = _Polynomial.from_bytes_wo_norm_and_wght(b=b)
        if len(b) == 8 + 4*decoded_poly.deg:
            return decoded_poly
        elif len(b) != 16 + 4*decoded_poly.deg:
            raise _DecodingError
        decoded_n_bytes: bytes = b[-8:-4]
        decoded_n: int = int.from_bytes(decoded_n_bytes, byteorder=ENDIANNESS)
        decoded_w_bytes: bytes = b[-4:]
        decoded_w: int = int.from_bytes(decoded_w_bytes, byteorder=ENDIANNESS)
        _, n, w = decoded_poly.coef_norm_wght
        if n != decoded_n or w != decoded_w:
            raise _DecodingError
        return decoded_poly
