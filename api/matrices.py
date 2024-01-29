from typing import List as _List, Tuple as _Tuple, TypeVar
from algebra.errors import _TYPE_MISMATCH_ERR
from algebra.matrices import _GeneralMatrix, _is_algebraic_class


MATRIX_STR_PREFIX: str = "Matrix"
T: type = TypeVar('T')


def is_algebraic_class(cls: type):
    return _is_algebraic_class(cls=cls)


class Matrix:
    matrix: _GeneralMatrix

    def __init__(self, matrix: _List[_List[T]]):
        self.matrix = _GeneralMatrix(matrix=matrix)

    @property
    def num_rows(self) -> int:
        return len(self.matrix)

    @property
    def num_cols(self) -> int:
        return len(self.matrix[0])

    @property
    def elem_class(self) -> type:
        return self.matrix.elem_class

    @property
    def coefs_norm_weight(self) -> _Tuple[_List[_List[_List[int]]], int, int]:
        return self.matrix.coefs_norm_weight

    def __str__(self):
        return MATRIX_STR_PREFIX + f"(matrix={self.matrix})"

    def __repr__(self):
        return self.__str__()

    def __len__(self):
        return len(self.matrix)

    # Don't expose the underlying _GeneralMatrix object
    def __iter__(self):
        return iter(self.matrix)

    def __getitem__(self, item):
        return self.matrix[item]

    def __setitem__(self, key, value):
        self.matrix[key] = value

    def __delitem__(self, key):
        self.matrix[key] = 0
    ##

    def __eq__(self, other):
        return isinstance(other, Matrix) and self.matrix == other.matrix

    def __add__(self, other):
        return Matrix(matrix=(self.matrix + other.matrix).matrix)

    def __neg__(self):
        return Matrix(matrix=(-self.matrix).matrix)

    def __sub__(self, other):
        return Matrix(matrix=(self.matrix - other.matrix).matrix)

    def __mul__(self, other):
        return Matrix(matrix=(self.matrix * other.matrix).matrix)

    ##
    def __mod__(self, other):
        return Matrix(matrix=self.matrix % other)
    ##