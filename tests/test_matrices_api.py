# test_matrices.py
import pytest
from algebra.matrices import GeneralMatrix as Matrix
from api.matrices import GeneralMatrix as ApiMatrix


@pytest.fixture
def sample_matrix():
    return ApiMatrix(Matrix([[1, 2], [3, 4]]))

def test_initialization_and_properties(sample_matrix):
    assert sample_matrix.rows == 2
    assert sample_matrix.cols == 2

def test_addition(sample_matrix):
    other = ApiMatrix(Matrix([[5, 6], [7, 8]]))
    result = sample_matrix + other
    expected = ApiMatrix(Matrix([[6, 8], [10, 12]]))
    assert len(result.vals.vals) == len(expected.vals.vals)
    for i in range(len(result.vals.vals)):
        assert len(result.vals.vals[i]) == len(expected.vals.vals[i])
        for j in range(len(result.vals.vals[i])):
            assert result.vals.vals[i][j] == expected.vals.vals[i][j]
    assert result.vals.vals == expected.vals.vals
    assert result.vals == expected.vals

def test_negation(sample_matrix):
    result = -sample_matrix
    expected = ApiMatrix(Matrix([[-1, -2], [-3, -4]]))
    assert result.vals == expected.vals

def test_subtraction(sample_matrix):
    other = ApiMatrix(Matrix([[5, 6], [7, 8]]))
    result = sample_matrix - other
    expected = ApiMatrix(Matrix([[-4, -4], [-4, -4]]))
    assert result.vals == expected.vals

def test_multiplication(sample_matrix):
    other = ApiMatrix(Matrix(((2, 0), (1, 2))))
    result = sample_matrix * other
    expected = ApiMatrix(Matrix(((4, 4), (10, 8))))
    assert result.vals == expected.vals

def test_string_representation(sample_matrix):
    assert str(sample_matrix) == "GeneralMatrix(matrix=GeneralMatrix(((1, 2), (3, 4))))"

def test_equality(sample_matrix):
    same_matrix = ApiMatrix(Matrix(((1, 2), (3, 4))))
    another_matrix = ApiMatrix(Matrix(((5, 6), (7, 8))))
    sample_is_same = sample_matrix == same_matrix
    assert sample_is_same
    assert sample_matrix != another_matrix
