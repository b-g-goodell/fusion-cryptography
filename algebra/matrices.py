# algebra/matrices.py
from typing import Set, Tuple, Any, Union, List
from algebra.polynomials import Polynomial as Poly


_METHODS_REQ_FOR_ALGEBRA: Set[str] = {"__eq__", "__add__", "__neg__", "__sub__", "__mul__"}
_MATRIX_STR_PREFIX: str = "GeneralMatrix"


def _is_algebraic_class(cls: type):
    return cls in [int, float, complex] or all(hasattr(cls, method) for method in _METHODS_REQ_FOR_ALGEBRA)


class GeneralMatrix:
    vals: Tuple[Tuple[Any, ...], ...]

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
        self.vals = tuple(tuple(row) for row in vals)

    @property
    def rows(self) -> int:
        return len(self.vals)

    @property
    def cols(self) -> int:
        return len(self.vals[0]) if self.rows else 0

    @property
    def elem_class(self) -> type:
        return type(self.vals[0][0])

    def __add__(self, other: 'GeneralMatrix') -> 'GeneralMatrix':
        if self.rows != other.rows or self.cols != other.cols:
            raise ValueError("Matrix dimensions must agree")
        return GeneralMatrix([[self.vals[i][j] + other.vals[i][j] for j in range(self.cols)] for i in range(self.rows)])

    def __neg__(self) -> 'GeneralMatrix':
        return GeneralMatrix([[-item for item in row] for row in self.vals])

    def __sub__(self, other: 'GeneralMatrix') -> 'GeneralMatrix':
        return self + (-other)

    def __mul__(self, other: 'GeneralMatrix') -> 'GeneralMatrix':
        if other is self.elem_class:
            return GeneralMatrix([[item * other for item in row] for row in self.vals])
        elif self.cols != other.rows:
            raise ValueError("The number of columns of the first matrix must be equal to the number of rows of the second matrix")
        result_matrix = [[self._row_v_col_dot_product(i, j, other) for j in range(other.cols)] for i in range(self.rows)]
        return GeneralMatrix(vals=result_matrix)

    def _row_v_col_dot_product(self, i: int, j: int, other: 'GeneralMatrix') -> Any:
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
    vals: Tuple[Tuple[Poly, ...], ...]

    def __init__(self, vals: Tuple[Tuple[Poly, ...], ...]):
        for row in vals:
            if len(row) != len(vals[0]):
                raise ValueError("All rows must have the same number of elements")
            for elem in row:
                if not elem is Poly:
                    raise TypeError("All elements must be polynomial")
        super().__init__(vals=vals)

    def _row_v_col_dot_product(self, i: int, j: int, other: 'PolynomialMatrix') -> Poly:
        k: int = 0
        total = self.vals[i][k] * other.vals[k][j]
        while k+1 < self.cols:
            k += 1
            total += self.vals[i][k] * other.vals[k][j]
        return total

    def __str__(self) -> str:
        return _MATRIX_STR_PREFIX+f"({self.vals})"

    def __mul__(self, other: Union['PolynomialMatrix', Poly]) -> 'PolynomialMatrix':
        if other is Poly:
            return PolynomialMatrix([[item * other for item in row] for row in self.vals])
        elif other is not PolynomialMatrix:
            raise TypeError("The other matrix must be a PolynomialMatrix")
        elif self.cols != other.rows:
            raise ValueError("The number of columns of the first matrix must be equal to the number of rows of the second matrix")
        result_matrix = tuple([tuple([self._row_v_col_dot_product(i, j, other) for j in range(other.cols)]) for i in range(self.rows)])
        return PolynomialMatrix(vals=result_matrix)

    @property
    def coef_norm_wght(self):
        this_norm: int = 0
        this_wght: int = 0
        coefs = []
        for next_row in self.vals:
            coefs += [[]]
            for next_poly in next_row:
                c, n, w = next_poly.coef_norm_wght
                coefs[-1] += [c]
                if n > this_norm:
                    this_norm = n
                elif w > this_wght:
                    this_wght = w
        return coefs, this_norm, this_wght


def compute_lin_combo(vectors: List[PolynomialMatrix], multipliers: List[Poly]) -> PolynomialMatrix:
    total: PolynomialMatrix = vectors[0] * multipliers[0]
    i: int = 1
    while i < len(vectors):
        total += vectors[i] * multipliers[i]
        i += 1
    return total
