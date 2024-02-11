import pytest
from algebra.errors import _DimensionMismatchError, _TypeMismatchError, _MustBeAlgebraicClassError
from algebra.ntt import _find_prou
from algebra.polynomials import _Polynomial as Poly
from algebra.matrices import _is_algebraic_class, _GeneralMatrix


MODULUS_FOR_TESTING: int = 17
DEGREE_FOR_TESTING: int = 8
ROOT_FOR_TESTING: int = _find_prou(mod=MODULUS_FOR_TESTING, deg=DEGREE_FOR_TESTING)
INV_ROOT_FOR_TESTING: int = pow(base=ROOT_FOR_TESTING, exp=MODULUS_FOR_TESTING - 2, mod=MODULUS_FOR_TESTING)
SAMPLE_SIZE: int = 2 ** 10


def test_is_algebraic_class():
    assert not _is_algebraic_class(cls="hello world".__class__())
    assert _is_algebraic_class(Poly)


# A mock class that supports the required arithmetic operations
class MockElement:
    def __init__(self, value):
        self.value = value

    def __add__(self, other):
        return MockElement(self.value + other.value)

    def __neg__(self):
        return MockElement(-self.value)

    def __sub__(self, other):
        return self.__add__(other.__neg__())

    def __mul__(self, other):
        return MockElement(self.value * other.value)

    def __eq__(self, other):
        if isinstance(other, MockElement):
            return self.value == other.value
        return False

    def __repr__(self):
        return f"MockElement({self.value})"


def test_general_matrix_addition():
    m1 = _GeneralMatrix(tuple([tuple([MockElement(1), MockElement(2)]), tuple([MockElement(3), MockElement(4)])]))
    m2 = _GeneralMatrix(tuple([tuple([MockElement(5), MockElement(6)]), tuple([MockElement(7), MockElement(8)])]))
    result = m1 + m2
    expected = _GeneralMatrix(tuple([tuple([MockElement(6), MockElement(8)]), tuple([MockElement(10), MockElement(12)])]))
    assert result._vals == expected._vals


def test_general_matrix_negation():
    m1 = _GeneralMatrix(tuple([tuple([MockElement(1), MockElement(-2)]), tuple([MockElement(-3), MockElement(4)])]))
    result = -m1
    expected = _GeneralMatrix(tuple([tuple([MockElement(-1), MockElement(2)]), tuple([MockElement(3), MockElement(-4)])]))
    assert result._vals == expected._vals


def test_general_matrix_subtraction():
    m1 = _GeneralMatrix(tuple([tuple([MockElement(10), MockElement(20)]), tuple([MockElement(30), MockElement(40)])]))
    m2 = _GeneralMatrix(tuple([tuple([MockElement(1), MockElement(2)]), tuple([MockElement(3), MockElement(4)])]))
    result = m1 - m2
    expected = _GeneralMatrix(tuple([tuple([MockElement(9), MockElement(18)]), tuple([MockElement(27), MockElement(36)])]))
    assert result._vals == expected._vals


def test_general_matrix_multiplication():
    m1 = _GeneralMatrix(tuple([tuple([MockElement(1), MockElement(2)]), tuple([MockElement(3), MockElement(4)])]))
    m2 = _GeneralMatrix(tuple([tuple([MockElement(5), MockElement(6)]), tuple([MockElement(7), MockElement(8)])]))
    result = m1 * m2
    expected = _GeneralMatrix(tuple([tuple([MockElement(19), MockElement(22)]), tuple([MockElement(43), MockElement(50)])]))
    assert result._vals == expected._vals


def test_matrix_dimension_mismatch_addition():
    m1 = _GeneralMatrix(tuple([tuple([MockElement(1)]), tuple([MockElement(2)])]))
    m2 = _GeneralMatrix(tuple([tuple([MockElement(1), MockElement(2)]), tuple([MockElement(3), MockElement(4)])]))
    with pytest.raises(_DimensionMismatchError):
        _ = m1 + m2


def test_matrix_dimension_mismatch_multiplication():
    m1 = _GeneralMatrix(tuple([tuple([MockElement(1), MockElement(2)])]))
    m2 = _GeneralMatrix(tuple([tuple([MockElement(3)]), tuple([MockElement(4)]), tuple([MockElement(5)])]))
    with pytest.raises(_DimensionMismatchError):
        _ = m1 * m2


def test_matrix_inconsistent_row_lengths():
    with pytest.raises(_DimensionMismatchError):
        _GeneralMatrix(tuple([tuple([MockElement(1), MockElement(2)]), tuple([MockElement(3)])]))


def test_matrix_different_element_types():
    with pytest.raises(_TypeMismatchError):
        _GeneralMatrix(tuple([tuple([MockElement(1), MockElement(2)]), (3, 4)]))


def test_matrix_element_type_lacking_operations():
    class BadMockElement:
        pass

    with pytest.raises(_MustBeAlgebraicClassError):
        _GeneralMatrix(tuple([tuple([BadMockElement(), BadMockElement()]), tuple([BadMockElement(), BadMockElement()])]))
