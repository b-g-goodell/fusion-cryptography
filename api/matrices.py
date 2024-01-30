from algebra.matrices import GeneralMatrix as Matrix, _MATRIX_STR_PREFIX


class GeneralMatrix(object):
    vals: Matrix

    def __init__(self, matrix: Matrix):
        self.vals = matrix

    @property
    def rows(self) -> int:
        return self.vals.rows

    @property
    def cols(self) -> int:
        return self.vals.cols

    @property
    def element_type(self) -> type:
        return self.vals.element_type

    def __eq__(self, other: 'GeneralMatrix') -> bool:
        return self.vals == other.vals

    def __add__(self, other: 'GeneralMatrix') -> 'GeneralMatrix':
        new_vals = self.vals + other.vals
        return GeneralMatrix(matrix=new_vals)

    def __neg__(self) -> 'GeneralMatrix':
        new_vals = -self.vals
        return GeneralMatrix(matrix=new_vals)

    def __sub__(self, other: 'GeneralMatrix') -> 'GeneralMatrix':
        return self + (-other)

    def __mul__(self, other: 'GeneralMatrix') -> 'GeneralMatrix':
        new_vals = self.vals * other.vals
        return GeneralMatrix(matrix=new_vals)

    def __str__(self) -> str:
        return _MATRIX_STR_PREFIX + f"(matrix={self.vals})"

    def __repr__(self) -> str:
        return self.__str__()
