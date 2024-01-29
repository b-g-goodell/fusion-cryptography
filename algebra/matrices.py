from copy import deepcopy
from typing import List as _List, Set as _Set, Tuple as _Tuple
from algebra.errors import _MUST_BE_INT_ERR, _MUST_BE_POS_INT_ERR, _MUST_BE_LIST_ERR, _MUST_BE_NONEMPTY_ERR, _DIMENSION_MISMATCH_ERR, _TYPE_MISMATCH_ERR, _COEFS_NORM_WEIGHT_METHOD_MISSING_ERR


_METHODS_REQ_FOR_ALGEBRA: _Set[str] = {"__eq__", "__add__", "__neg__", "__sub__", "__mul__"}
_MATRIX_STR_PREFIX: str = "_GeneralMatrix"


def _is_algebraic_class(cls: type):
    return all(hasattr(cls, method) for method in _METHODS_REQ_FOR_ALGEBRA)


class _GeneralMatrix:
    elem_class: type
    matrix: _List[_List[object]]

    def __init__(self, matrix: _List[list]):
        if not isinstance(matrix, list) or any(not isinstance(x, list) for x in matrix):
            raise ValueError(_MUST_BE_LIST_ERR)
        elif not matrix or any(not x for x in matrix):
            raise ValueError(_MUST_BE_NONEMPTY_ERR)
        elif not all(len(row) == len(matrix[0]) for row in matrix):
            raise ValueError(_DIMENSION_MISMATCH_ERR)
        elif not _is_algebraic_class(cls=matrix[0][0].__class__):
            raise ValueError(_TYPE_MISMATCH_ERR)
        elif not all(
            all(isinstance(item, matrix[0][0].__class__) for item in row)
            for row in matrix
        ):
            raise ValueError(_TYPE_MISMATCH_ERR)
        self.elem_class = matrix[0][0].__class__
        self.matrix = matrix

    def __str__(self):
        return _MATRIX_STR_PREFIX + f"(elem_class={self.elem_class}, matrix={self.matrix})"

    def __repr__(self):
        return self.__str__()

    ##
    def __len__(self):
        return len(self.matrix)

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
        return isinstance(other, _GeneralMatrix) and self.elem_class == other.elem_class and self.matrix == other.matrix

    def __add__(self, other):
        if other == 0:
            return self
        elif not isinstance(other, _GeneralMatrix) or self.elem_class != other.elem_class:
            raise TypeError(_TYPE_MISMATCH_ERR)
        elif len(self.matrix) != len(other.matrix) or len(self.matrix[0]) != len(other.matrix[0]):
            raise TypeError(_DIMENSION_MISMATCH_ERR)
        resulting_matrix = [
            [self.matrix[i][j] + other.matrix[i][j] for j in range(len(self.matrix[0]))]
            for i in range(len(self.matrix))
        ]
        return _GeneralMatrix(matrix=resulting_matrix)

    def __radd__(self, other):
        if other == 0:
            return self
        return self + other

    def __neg__(self):
        resulting_matrix = [
            [-self.matrix[i][j] for j in range(len(self.matrix[0]))]
            for i in range(len(self.matrix))
        ]
        return _GeneralMatrix(matrix=resulting_matrix)

    def __sub__(self, other):
        return self + (-other)

    def __mul__(self, other):
        if other == 0:
            return 0
        elif other == 1:
            return self
        elif not isinstance(other, self.elem_class) and not (isinstance(other, _GeneralMatrix) and self.elem_class == other.elem_class):
            raise TypeError(_TYPE_MISMATCH_ERR)
        elif isinstance(other, self.elem_class):
            # scalar multiplication
            resulting_matrix = [
                [self.matrix[i][j] * other for j in range(len(self.matrix[0]))]
                for i in range(len(self.matrix))
            ]
            return _GeneralMatrix(matrix=resulting_matrix)
        elif len(self.matrix) != len(other.matrix) or len(self.matrix[0]) != len(other.matrix[0]):
            raise TypeError(_DIMENSION_MISMATCH_ERR)
        result: _GeneralMatrix = deepcopy(self)
        result.matrix = [
            [0 for j in range(len(other.matrix[0]))] for i in range(len(self.matrix))
        ]
        for i in range(len(self.matrix)):
            for j in range(len(other.matrix[0])):
                next_data = self.matrix[i][0] * other.matrix[0][j]
                for k in range(1, len(self.matrix[0])):
                    next_data += self.matrix[i][k] * other.matrix[k][j]
                result.matrix[i][j] = next_data
        return result

    def __mod__(self, other):
        if not isinstance(other, int):
            raise TypeError(_MUST_BE_INT_ERR)
        elif other <= 1:
            raise ValueError(_MUST_BE_POS_INT_ERR)
        elem_class: type = self.elem_class
        resulting_matrix: _List[_List[elem_class]] = [
            [self.matrix[i][j] % other for j in range(len(self.matrix[0]))]
            for i in range(len(self.matrix))
        ]
        return _GeneralMatrix(matrix=resulting_matrix)

    @property
    def coefs_norm_weight(self) -> _Tuple[_List[_List[_List[int]]], int, int]:
        if not all(hasattr(z, "coefs_norm_weight") for y in self.matrix for z in y):
            raise NotImplementedError(_COEFS_NORM_WEIGHT_METHOD_MISSING_ERR)
        cnw: _List[_List[_Tuple[_List[int], int, int]]] = [[z.coefs_norm_weight for z in y] for y in self.matrix]
        coefs = [[z[0] for z in y] for y in cnw]
        norm = max(max(z[1] for z in y) for y in cnw)
        weight = max(max(z[2] for z in y) for y in cnw)
        return coefs, norm, weight
