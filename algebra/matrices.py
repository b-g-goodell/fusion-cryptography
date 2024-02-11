"""
algebra/matrices.py

Handle linear algebra of matrices whose elements inside support algebraic functions like addition, negation, sub-
traction, and multiplication. In particular, handle  matrices of polynomials from Z_p[X]/(X^d+1) for p a prime and d a
power of 2 such that (p-1)%(2d)==0.
"""
# TODO: add ser/deser
from algebra.errors import _DimensionMismatchError, _TypeMismatchError, _MustBeAlgebraicClassError, _MustBeIntError
from algebra.polynomials import _Polynomial as _Poly, ENDIANNESS, _DecodingError

_METHODS_REQ_FOR_ALGEBRA: set[str] = {"__eq__", "__add__", "__neg__", "__sub__", "__mul__"}
_MATRIX_STR_PREFIX: str = "GeneralMatrix"


def _is_algebraic_class(cls: type):
    return (cls in [int, float, complex] or
            all(hasattr(cls, method) for method in _METHODS_REQ_FOR_ALGEBRA))


class _GeneralMatrix:
    """
    A class to represent a general matrix whose elements support algebraic operations like addition, negation,
    subtraction, and multiplication. The elements must be of the same type and must support the required operations.

    :ivar _vals: The tuple of tuples representing the matrix
    """
    _vals: tuple[tuple[..., ...], ...]

    def __init__(self, vals: tuple[tuple[..., ...], ...]):
        # Check for consistent row lengths and element types
        element_type: type = type(None)
        for row in vals:
            if len(row) != len(vals[0]):
                raise _DimensionMismatchError
            for elem in row:
                if element_type is type(None):
                    element_type = type(elem)
                elif not (type(elem) is element_type):
                    raise _TypeMismatchError
        # Checking for required operation support in the element type
        if not all(hasattr(element_type, op) for op in _METHODS_REQ_FOR_ALGEBRA):
            raise _MustBeAlgebraicClassError
        self._vals = tuple(tuple(row) for row in vals)

    @property
    def rows(self) -> int:
        """ Return the number of rows """
        return len(self._vals)

    @property
    def cols(self) -> int:
        """ Return the number of columns """
        return len(self._vals[0]) if self.rows else 0

    @property
    def elem_class(self) -> type:
        """ Return the class of the elements """
        return type(self._vals[0][0])

    @property
    def vals(self) -> tuple[tuple[..., ...], ...]:
        """ Return the matrix as a tuple of tuples """
        return self._vals

    def __add__(self, other: '_GeneralMatrix') -> '_GeneralMatrix':
        if self.rows != other.rows or self.cols != other.cols:
            raise _DimensionMismatchError
        return _GeneralMatrix(tuple([tuple([self.vals[i][j] + other.vals[i][j] for j in range(self.cols)]) for i in range(self.rows)]))

    def __neg__(self) -> '_GeneralMatrix':
        return _GeneralMatrix(tuple([tuple([-item for item in row]) for row in self.vals]))

    def __sub__(self, other: '_GeneralMatrix') -> '_GeneralMatrix':
        return self + (-other)

    def __mul__(self, other: '_GeneralMatrix') -> '_GeneralMatrix':
        if other is self.elem_class:
            return _GeneralMatrix(tuple([tuple([item * other for item in row]) for row in self.vals]))
        elif self.cols != other.rows:
            raise _DimensionMismatchError
        result_matrix = tuple([tuple([self._row_v_col_dot_product(i, j, other) for j in range(other.cols)]) for i in range(self.rows)])
        return _GeneralMatrix(vals=result_matrix)

    def _row_v_col_dot_product(self, i: int, j: int, other: '_GeneralMatrix'):
        """ Returns the i-th row of self dotted with the j-th column of other. """
        if not isinstance(other, _GeneralMatrix):
            raise _TypeMismatchError
        k: int = 0
        total = self.vals[i][k] * other.vals[k][j]
        while k+1 < self.cols:
            k += 1
            total += self.vals[i][k] * other.vals[k][j]
        return total

    def __str__(self) -> str:
        return _MATRIX_STR_PREFIX+f"({self.vals})"

    def __mod__(self, other):
        if not isinstance(other, int):
            raise _MustBeIntError
        return _GeneralMatrix(tuple([tuple([item % other for item in row]) for row in self.vals]))

    def __eq__(self, other: '_GeneralMatrix') -> bool:
        return self.vals == other.vals


class _PolynomialMatrix(_GeneralMatrix):
    """
    A class to represent a matrix of polynomials from Z_p[X]/(X^d+1).

    :ivar _vals: The tuple of tuples representing the matrix
    :ivar _coefs: The tuple of tuples representing the coefficients of the polynomials
    :ivar _norm: The maximum absolute value of the coefficients
    :ivar _wght: The number of non-zero coefficients
    """
    _vals: tuple[tuple[_Poly, ...], ...]

    def __init__(self, vals: tuple[tuple[_Poly, ...], ...]):
        for row in vals:
            if len(row) != len(vals[0]):
                raise _DimensionMismatchError
            for elem in row:
                if not isinstance(elem, _Poly):
                    raise _TypeMismatchError
        super().__init__(vals=vals)

    def _row_v_col_dot_product(self, i: int, j: int, other: '_PolynomialMatrix') -> _Poly:
        if not isinstance(other, _PolynomialMatrix):
            raise _TypeMismatchError
        return super()._row_v_col_dot_product(i, j, other)

    def __str__(self) -> str:
        return _MATRIX_STR_PREFIX+f"({self.vals})"

    def __mul__(self, other) -> '_PolynomialMatrix':
        if isinstance(other, _Poly):
            return _PolynomialMatrix(tuple([tuple([item * other for item in row]) for row in self.vals]))
        elif not isinstance(other, _PolynomialMatrix):
            raise _TypeMismatchError
        elif self.cols != other.rows:
            raise _DimensionMismatchError
        result_matrix = tuple([tuple([
            self._row_v_col_dot_product(i, j, other) for j in range(other.cols)]) for i in range(self.rows)])
        return _PolynomialMatrix(vals=result_matrix)

    @property
    def coef_norm_wght(self) -> tuple[tuple[tuple[tuple[int, ...], ...], ...], int, int]:
        this_norm: int = 0
        this_wght: int = 0
        coefs_list: list[list[tuple[int, ...]]] = []
        c: tuple[int, ...]
        n: int
        w: int
        for next_row in self._vals:
            coefs_list += [[]]
            for next_elem in next_row:
                c, n, w = next_elem.coef_norm_wght
                coefs_list[-1] += [c]
                if n > this_norm:
                    this_norm = n
                elif w > this_wght:
                    this_wght = w
        return tuple([tuple([next_col for next_col in next_row]) for next_row in coefs_list]), this_norm, this_wght

    def to_bytes(self) -> bytes:
        q: int = self.vals[0][0].modulus
        d: int = self.vals[0][0].degree
        num_rows: int = self.rows
        num_cols: int = self.cols
        x: bytes = num_rows.to_bytes(length=4, byteorder=ENDIANNESS) + num_cols.to_bytes(length=4, byteorder=ENDIANNESS)
        x += q.to_bytes(length=4, byteorder=ENDIANNESS) + d.to_bytes(length=4, byteorder=ENDIANNESS)
        for row in self.vals:
            for next_poly in row:
                next_poly_bytes: bytes = next_poly.to_bytes()
                x += next_poly_bytes[8:8+4*d]
        _, n, w = self.coef_norm_wght
        x += n.to_bytes(length=4, byteorder=ENDIANNESS) + w.to_bytes(length=4, byteorder=ENDIANNESS)
        return x

    @staticmethod
    def from_bytes_wo_norm_and_wght(b: bytes):
        if len(b) < 16:
            raise _DecodingError
        num_rows_bytes: bytes = b[:4]
        num_rows: int = int.from_bytes(num_rows_bytes, byteorder=ENDIANNESS)
        num_cols_bytes: bytes = b[4:8]
        num_cols: int = int.from_bytes(num_cols_bytes, byteorder=ENDIANNESS)
        mod_bytes: bytes = b[8:12]
        modulus: int = int.from_bytes(mod_bytes, byteorder=ENDIANNESS)
        deg_bytes: bytes = b[12:16]
        degree: int = int.from_bytes(deg_bytes, byteorder=ENDIANNESS)
        if len(b) < 16+4*degree*num_rows*num_cols:
            raise _DecodingError
        vals_bytes: bytes = b[16:16+4*degree*num_rows*num_cols]
        vals: list[list[_Poly]] = []
        for i in range(num_rows):
            vals += [[]]
            for j in range(num_cols):
                next_poly_bytes: bytes = vals_bytes[4*degree*(i*num_cols+j):4*degree*(i*num_cols+j)+4*degree]
                vals[-1] += [_Poly.from_bytes_wo_mod_and_deg(b=next_poly_bytes)]
        return _PolynomialMatrix(vals=tuple([tuple(row) for row in vals]))

    @staticmethod
    def from_bytes(b: bytes):
        decoded_matrix: _PolynomialMatrix = _PolynomialMatrix.from_bytes_wo_norm_and_wght(b=b)
        num_rows: int = decoded_matrix.rows
        num_cols: int = decoded_matrix.cols
        degree: int = decoded_matrix.vals[0][0].degree
        if len(b) == 16+4*degree*num_rows*num_cols:
            return decoded_matrix
        elif len(b) != 24+4*degree*num_rows*num_cols:
            raise _DecodingError
        norm_bytes: bytes = b[16+4*degree*num_rows*num_cols:20+4*degree*num_rows*num_cols]
        decoded_norm: int = int.from_bytes(norm_bytes, byteorder=ENDIANNESS)
        wght_bytes: bytes = b[20+4*degree*num_rows*num_cols+4:24+4*degree*num_rows*num_cols]
        decoded_wght: int = int.from_bytes(wght_bytes, byteorder=ENDIANNESS)
        _, n, w = decoded_matrix.coef_norm_wght
        if n != decoded_norm or w != decoded_wght:
            raise _DecodingError
        return decoded_matrix


def _compute_lin_combo(vectors: list[_PolynomialMatrix], multipliers: list[_Poly]) -> _PolynomialMatrix:
    if len(vectors) != len(multipliers):
        raise _DimensionMismatchError
    n: int = len(vectors)
    total: _PolynomialMatrix = vectors[0] * multipliers[0]
    i: int = 1
    while i < n:
        total += vectors[i] * multipliers[i]
        i += 1
    return total
