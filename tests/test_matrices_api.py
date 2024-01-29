import pytest
from random import randrange
from typing import List
from api.matrices import Matrix, is_algebraic_class
from api.ntt import find_prou
from api.polynomials import Polynomial as Poly
from copy import deepcopy


MODULUS_FOR_TESTING: int = 17
DEGREE_FOR_TESTING: int = 8
ROOT_FOR_TESTING: int = find_prou(mod=MODULUS_FOR_TESTING, deg=DEGREE_FOR_TESTING)
INV_ROOT_FOR_TESTING: int = pow(base=ROOT_FOR_TESTING, exp=MODULUS_FOR_TESTING - 2, mod=MODULUS_FOR_TESTING)
SAMPLE_SIZE: int = 2 ** 10


def test_is_algebraic_class():
    assert not is_algebraic_class(cls="hello world".__class__())
    assert is_algebraic_class(Poly)


def test_general_matrix():
    for _ in range(SAMPLE_SIZE):
        # Generate random left-matrix [[a_left, b_left], [c_left, d_left]]
        # and right-matrix [[a_right, b_right], [c_right, d_right]]
        a_left_coef: int = randrange(1, MODULUS_FOR_TESTING)
        a_right_coef: int = randrange(1, MODULUS_FOR_TESTING)
        a_left_index: int = randrange(DEGREE_FOR_TESTING)
        a_right_index: int = randrange(DEGREE_FOR_TESTING)
        a_left_coefs: List[int] = [
            0 if i != a_left_index else a_left_coef for i in range(DEGREE_FOR_TESTING)
        ]
        a_right_coefs: List[int] = [
            0 if i != a_right_index else a_right_coef for i in range(DEGREE_FOR_TESTING)
        ]
        a_left: Poly = Poly(
            modulus=MODULUS_FOR_TESTING,
            values=a_left_coefs,
            representation="coefficient",
        )
        a_right: Poly = Poly(
            modulus=MODULUS_FOR_TESTING,
            values=a_right_coefs,
            representation="coefficient",
        )

        b_left_coef: int = randrange(1, MODULUS_FOR_TESTING)
        b_right_coef: int = randrange(1, MODULUS_FOR_TESTING)
        b_left_index: int = randrange(DEGREE_FOR_TESTING)
        b_right_index: int = randrange(DEGREE_FOR_TESTING)
        b_left_coefs: List[int] = [
            0 if i != b_left_index else b_left_coef for i in range(DEGREE_FOR_TESTING)
        ]
        b_right_coefs: List[int] = [
            0 if i != b_right_index else b_right_coef for i in range(DEGREE_FOR_TESTING)
        ]
        b_left: Poly = Poly(
            modulus=MODULUS_FOR_TESTING,
            values=b_left_coefs,
            representation="coefficient",
        )
        b_right: Poly = Poly(
            modulus=MODULUS_FOR_TESTING,
            values=b_right_coefs,
            representation="coefficient",
        )

        c_left_coef: int = randrange(1, MODULUS_FOR_TESTING)
        c_right_coef: int = randrange(1, MODULUS_FOR_TESTING)
        c_left_index: int = randrange(DEGREE_FOR_TESTING)
        c_right_index: int = randrange(DEGREE_FOR_TESTING)
        c_left_coefs: List[int] = [
            0 if i != c_left_index else c_left_coef for i in range(DEGREE_FOR_TESTING)
        ]
        c_right_coefs: List[int] = [
            0 if i != c_right_index else c_right_coef for i in range(DEGREE_FOR_TESTING)
        ]

        c_left: Poly = Poly(
            modulus=MODULUS_FOR_TESTING,
            values=c_left_coefs,
            representation="coefficient",
        )
        c_right: Poly = Poly(
            modulus=MODULUS_FOR_TESTING,
            values=c_right_coefs,
            representation="coefficient",
        )

        d_left_coef: int = randrange(1, MODULUS_FOR_TESTING)
        d_right_coef: int = randrange(1, MODULUS_FOR_TESTING)
        d_left_index: int = randrange(DEGREE_FOR_TESTING)
        d_right_index: int = randrange(DEGREE_FOR_TESTING)
        d_left_coefs: List[int] = [
            0 if i != d_left_index else d_left_coef for i in range(DEGREE_FOR_TESTING)
        ]
        d_right_coefs: List[int] = [
            0 if i != d_right_index else d_right_coef for i in range(DEGREE_FOR_TESTING)
        ]
        d_left: Poly = Poly(
            modulus=MODULUS_FOR_TESTING,
            values=d_left_coefs,
            representation="coefficient",
        )
        d_right: Poly = Poly(
            modulus=MODULUS_FOR_TESTING,
            values=d_right_coefs,
            representation="coefficient",
        )

        left_matrix: Matrix = Matrix(
            matrix=[
                [deepcopy(a_left), deepcopy(b_left)],
                [deepcopy(c_left), deepcopy(d_left)],
            ]
        )
        right_matrix: Matrix = Matrix(
            matrix=[
                [deepcopy(a_right), deepcopy(b_right)],
                [deepcopy(c_right), deepcopy(d_right)],
            ]
        )

        # Test the left-matrix
        assert left_matrix.matrix[0][0] == a_left
        assert left_matrix.matrix[0][1] == b_left
        assert left_matrix.matrix[1][0] == c_left
        assert left_matrix.matrix[1][1] == d_left
        assert left_matrix.elem_class == Poly

        # Test the right-matrix
        assert right_matrix.matrix[0][0] == a_right
        assert right_matrix.matrix[0][1] == b_right
        assert right_matrix.matrix[1][0] == c_right
        assert right_matrix.matrix[1][1] == d_right
        assert right_matrix.elem_class == Poly

        # Test left-matrix + right-matrix
        expected_sum: Matrix = Matrix(
            matrix=[
                [a_left + a_right, b_left + b_right],
                [c_left + c_right, d_left + d_right],
            ]
        )
        observed_sum: Matrix = left_matrix + right_matrix
        assert observed_sum == expected_sum

        # Test the left-matrix * right-matrix
        expected_product: Matrix = Matrix(
            matrix=[
                [
                    a_left * a_right + b_left * c_right,
                    a_left * b_right + b_left * d_right,
                ],
                [
                    c_left * a_right + d_left * c_right,
                    c_left * b_right + d_left * d_right,
                ],
            ]
        )
        observed_product: Matrix = left_matrix * right_matrix
        assert observed_product == expected_product
