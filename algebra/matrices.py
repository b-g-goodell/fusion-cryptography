# algebra/matrices.py
from typing import Set, Tuple, Any, Union, List
from algebra.polynomials import Polynomial as Poly


_METHODS_REQ_FOR_ALGEBRA: Set[str] = {"__eq__", "__add__", "__neg__", "__sub__", "__mul__"}
_MATRIX_STR_PREFIX: str = "GeneralMatrix"


def _is_algebraic_class(cls: type):
    return cls in [int, float, complex] or all(hasattr(cls, method) for method in _METHODS_REQ_FOR_ALGEBRA)


class GeneralMatrix:
    _vals: Tuple[Tuple[Any, ...], ...]

    def __init__(self, vals: Tuple[Tuple[Any, ...], ...]):
        # Check for consistent row lengths and element types
        element_type: type = type(None)
        for row in vals:
            if len(row) != len(vals[0]):
                raise ValueError("All rows must have the same number of elements")
            for elem in row:
                if element_type is type(None):
                    element_type = type(elem)
                elif not (type(elem) is element_type):
                    raise TypeError("All elements must be of the same type")
        # Checking for required operation support in the element type
        required_ops = ('__add__', '__neg__', '__sub__', '__mul__')
        if not all(hasattr(element_type, op) for op in required_ops):
            raise TypeError("Element type must support add, neg, sub, and mul operations")
        self._vals = tuple(tuple(row) for row in vals)

    @property
    def rows(self) -> int:
        return len(self._vals)

    @property
    def cols(self) -> int:
        return len(self._vals[0]) if self.rows else 0

    @property
    def elem_class(self) -> type:
        return type(self._vals[0][0])

    @property
    def vals(self) -> Tuple[Tuple[Any, ...], ...]:
        return self._vals

    def __add__(self, other: 'GeneralMatrix') -> 'GeneralMatrix':
        if self.rows != other.rows or self.cols != other.cols:
            raise ValueError("Matrix dimensions must agree")
        return GeneralMatrix(tuple([tuple([self.vals[i][j] + other.vals[i][j] for j in range(self.cols)]) for i in range(self.rows)]))

    def __neg__(self) -> 'GeneralMatrix':
        return GeneralMatrix(tuple([tuple([-item for item in row]) for row in self.vals]))

    def __sub__(self, other: 'GeneralMatrix') -> 'GeneralMatrix':
        return self + (-other)

    def __mul__(self, other: 'GeneralMatrix') -> 'GeneralMatrix':
        if other is self.elem_class:
            return GeneralMatrix(tuple([tuple([item * other for item in row]) for row in self.vals]))
        elif self.cols != other.rows:
            raise ValueError("The number of columns of the first matrix must be equal to the number of rows of the second matrix")
        result_matrix = tuple([tuple([self._row_v_col_dot_product(i, j, other) for j in range(other.cols)]) for i in range(self.rows)])
        return GeneralMatrix(vals=result_matrix)

    def _row_v_col_dot_product(self, i: int, j: int, other: 'GeneralMatrix') -> Any:
        if not isinstance(other, GeneralMatrix):
            raise TypeError("The other matrix must be a GeneralMatrix")
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
            raise TypeError("The modulus must be an integer")
        return GeneralMatrix(tuple([tuple([item % other for item in row]) for row in self.vals]))

    def __eq__(self, other: 'GeneralMatrix') -> bool:
        return self.vals == other.vals


class PolynomialMatrix(GeneralMatrix):
    _vals: Tuple[Tuple[Poly, ...], ...]
    _coefs: Tuple[Tuple[Tuple[int, ...], ...], ...]
    _norm: int
    _wght: int

    def __init__(self, vals: Tuple[Tuple[Poly, ...], ...]):
        for row in vals:
            if len(row) != len(vals[0]):
                raise ValueError("All rows must have the same number of elements")
            for elem in row:
                if not elem is Poly:
                    raise TypeError("All elements must be polynomial")
        super().__init__(vals=vals)

    def _row_v_col_dot_product(self, i: int, j: int, other: 'PolynomialMatrix') -> Poly:
        if not isinstance(other, PolynomialMatrix):
            raise TypeError("The other matrix must be a PolynomialMatrix")
        return super()._row_v_col_dot_product(i, j, other)

    def __str__(self) -> str:
        return _MATRIX_STR_PREFIX+f"({self.vals})"

    def __mul__(self, other: Union['PolynomialMatrix', Poly]) -> 'PolynomialMatrix':
        if other is Poly:
            return PolynomialMatrix(tuple([tuple([item * other for item in row]) for row in self.vals]))
        elif other is not PolynomialMatrix:
            raise TypeError("The other matrix must be a PolynomialMatrix")
        elif self.cols != other.rows:
            raise ValueError("The number of columns of the first matrix must be equal to the number of rows of the second matrix")
        result_matrix = tuple([tuple([self._row_v_col_dot_product(i, j, other) for j in range(other.cols)]) for i in range(self.rows)])
        return PolynomialMatrix(vals=result_matrix)

    @property
    def coef_norm_wght(self) -> Tuple[Tuple[Tuple[Tuple[int, ...], ...], ...], int, int]:
        this_norm: int = 0
        this_wght: int = 0
        coefs_list: List[List[Tuple[int, ...]]] = []
        c: Tuple[int, ...]
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
        self._coefs = tuple([tuple(row) for row in coefs_list])
        self._norm = this_norm
        self._wght = this_wght
        return self._coefs, self._norm, self._wght

    @property
    def coefs(self) -> Tuple[Tuple[Tuple[int, ...], ...], ...]:
        return self._coefs

    @property
    def norm(self) -> int:
        return self._norm

    @property
    def wght(self) -> int:
        return self._wght


def compute_lin_combo(vectors: List[PolynomialMatrix], multipliers: List[Poly]) -> PolynomialMatrix:
    if len(vectors) != len(multipliers):
        raise ValueError("The number of vectors and multipliers must be the same")
    n: int = len(vectors)
    total: PolynomialMatrix = vectors[0] * multipliers[0]
    i: int = 1
    while i < n:
        total += vectors[i] * multipliers[i]
        i += 1
    return total
