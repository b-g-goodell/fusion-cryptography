import pytest
from api.polynomials import Polynomial
from api.errors import _INVALID_REP_TYPE_ERR

# Constants for testing
MODULUS = 17
COEFFICIENTS = tuple(list(range(1,9)))
COEFFICIENTS_REP = 'coefficient'
NTT_REP = 'ntt'
INVALID_REP = 'invalid_rep'
ZERO_POLYNOMIAL_VALUES = tuple([0] * len(COEFFICIENTS))

# Mocks for external dependencies not provided
_POLYNOMIAL_REPRESENTATION_TYPES = [COEFFICIENTS_REP, NTT_REP]


@pytest.fixture
def poly_coefficients():
    return Polynomial(MODULUS, COEFFICIENTS, COEFFICIENTS_REP)


@pytest.fixture
def poly_ntt():
    return Polynomial(MODULUS, COEFFICIENTS, NTT_REP)


@pytest.fixture
def zero_polynomial():
    return Polynomial(MODULUS, ZERO_POLYNOMIAL_VALUES, NTT_REP)


def test_initialization_with_coefficients():
    p = Polynomial(MODULUS, COEFFICIENTS, COEFFICIENTS_REP)
    c,n,w = p.coefs_norm_weight
    assert isinstance(p, Polynomial)
    assert p.mod == MODULUS
    assert c == COEFFICIENTS  # Assuming transformation to NTT keeps values same in this example.


def test_initialization_with_ntt():
    p = Polynomial(MODULUS, COEFFICIENTS, NTT_REP)
    assert isinstance(p, Polynomial)
    assert p.mod == MODULUS
    assert p.vals == COEFFICIENTS


def test_initialization_with_invalid_representation():
    with pytest.raises(ValueError) as e:
        Polynomial(MODULUS, COEFFICIENTS, INVALID_REP)
    assert str(e.value) == _INVALID_REP_TYPE_ERR


def test_polynomial_properties(poly_ntt):
    assert poly_ntt.mod == MODULUS
    # Further properties can be tested similar to above, if relevant getters are implemented.


def test_polynomial_equality(poly_coefficients, poly_ntt):
    p = Polynomial(MODULUS, COEFFICIENTS, COEFFICIENTS_REP)
    assert poly_coefficients == p
    assert poly_coefficients != poly_ntt  # Assume different representation means different polynomials.


def test_polynomial_addition(poly_coefficients, poly_ntt, zero_polynomial):
    result = poly_coefficients + poly_ntt
    assert isinstance(result, Polynomial)
    # The result should be tested further based on the known behavior of addition operation


def test_polynomial_subtraction(poly_coefficients, zero_polynomial):
    result = poly_coefficients - zero_polynomial
    assert isinstance(result, Polynomial)
    # Subtraction result tests similar to addition


def test_polynomial_multiplication(poly_coefficients, zero_polynomial):
    result = poly_coefficients * zero_polynomial
    assert isinstance(result, Polynomial)
    # Multiplication result tests similar to addition and subtraction


def test_polynomial_negation(poly_coefficients):
    negated = -poly_coefficients
    assert isinstance(negated, Polynomial)
    # The negation result should be tested based on the behavior of negation operation


def test_string_representation(poly_ntt):
    expected_representation = f"Polynomial(ntt={poly_ntt.ntt_rep})"
    assert str(poly_ntt) == expected_representation


def test_representation_equality():
    p1 = Polynomial(MODULUS, COEFFICIENTS, COEFFICIENTS_REP)
    p2 = Polynomial(MODULUS, COEFFICIENTS, COEFFICIENTS_REP)
    assert repr(p1) == repr(p2)