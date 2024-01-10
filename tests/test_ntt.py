import pytest
from typing import List, Tuple
from secrets import randbelow
from algebra.ntt import (
    is_odd_prime,
    is_pow_two_geq_two,
    has_primitive_root_of_unity,
    bit_reverse_copy,
    cent,
    is_root_of_unity,
    is_primitive_root,
    find_primitive_root,
    cooley_tukey_ntt,
    gentleman_sande_intt,
    ntt_poly_mult,
)

SAMPLE_SIZE: int = 2**5
IS_ODD_PRIME_TEST_DATA = [
    (3, True),
    (4, False),
    (5, True),
    (6, False),
    (29, True),
    (30, False),
    (1, False),
    (2, False),
]
IS_POW_TWO_GEQ_TWO_TEST_DATA = [
    (3, False),
    (2, True),
    (1, False),
    (4, True),
    (8, True),
    (7, False)
]
LOG2_D_MIN: int = 2
LOG2_D_MAX: int = 6
Q_MAX: int = 2**17
PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE: List[Tuple[int, int]] = []
for log2d in range(LOG2_D_MIN, LOG2_D_MAX + 1):
    tmp_d: int = 1 << log2d
    tmp_q: int = 2 * tmp_d + 1
    while tmp_q < Q_MAX:
        while not is_odd_prime(tmp_q) and tmp_q < Q_MAX:
            tmp_q += 2 * tmp_d
        if is_odd_prime(tmp_q) and tmp_q < Q_MAX:
            PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE.append((tmp_d, tmp_q))
            find_primitive_root(modulus=tmp_q, root_order=2 * tmp_d)
            tmp_q *= 2
            tmp_q -= (tmp_q - 1) % (2 * tmp_d)
            assert (tmp_q - 1) % (2 * tmp_d) == 0
HAS_PRIMITIVE_ROOT_OF_UNITY_TEST_DATA = [
    (17, 2, True),
    (17, 4, True),
    (17, 8, True),
    (17, 16, True),
    (17, 3, False),
    (17, 5, False),
    (17, 6, False),
    (17, 7, False)
]
BIT_REVERSE_COPY_TEST_DATA = [
    ([0, 1, 2, 3], [0, 2, 1, 3]),
    ([0, 1, 2, 3, 4, 5, 6, 7], [0, 4, 2, 6, 1, 5, 3, 7])
]
CENT_TEST_DATA = [
    (0, 17, 8, 5, 0),
    (1, 17, 8, 5, 1),
    (2, 17, 8, 5, 2),
    (3, 17, 8, 5, 3),
    (4, 17, 8, 5, 4),
    (5, 17, 8, 5, 5),
    (6, 17, 8, 5, 6),
    (7, 17, 8, 5, 7),
    (8, 17, 8, 5, 8),
    (9, 17, 8, 5, -8),
    (10, 17, 8, 5, -7),
    (11, 17, 8, 5, -6),
    (12, 17, 8, 5, -5),
    (13, 17, 8, 5, -4),
    (14, 17, 8, 5, -3),
    (15, 17, 8, 5, -2),
    (16, 17, 8, 5, -1),
]
IS_ROOT_OF_UNITY_TEST_DATA = [
    (3, 17, 16, True),
    (2, 17, 16, True),
]
IS_PRIMITIVE_ROOT_OF_UNITY_TEST_DATA = [
    (3, 17, 16, True),
    (2, 17, 16, False),
]
FIND_PRIMITIVE_ROOT_TEST_DATA = [
    (17, 16, 3),
    (257, 128, 9)
]
COOLEY_TUKEY_NTT_TEST_DATA = [
    (17, 8, 5, 16, find_primitive_root(modulus=17, root_order=16),
     pow(find_primitive_root(modulus=17, root_order=16), 15, 17),
     [pow(find_primitive_root(modulus=17, root_order=16), i, 17) for i in range(16)],
     bit_reverse_copy([pow(find_primitive_root(modulus=17, root_order=16), i, 17) for i in range(16)]),
     [1, 0, 0, 0, 0, 0, 0, 0], [1, 1, 1, 1, 1, 1, 1, 1]),
    (17, 8, 5, 16, find_primitive_root(modulus=17, root_order=16),
     pow(find_primitive_root(modulus=17, root_order=16), 15, 17),
     [pow(find_primitive_root(modulus=17, root_order=16), i, 17) for i in range(16)],
     bit_reverse_copy([pow(find_primitive_root(modulus=17, root_order=16), i, 17) for i in range(16)]),
     [0, 1, 0, 0, 0, 0, 0, 0], [-8, 8, 8, -8, -2, 2, 2, -2]),
    (17, 8, 5, 16, find_primitive_root(modulus=17, root_order=16),
     pow(find_primitive_root(modulus=17, root_order=16), 15, 17),
     [pow(find_primitive_root(modulus=17, root_order=16), i, 17) for i in range(16)],
     bit_reverse_copy([pow(find_primitive_root(modulus=17, root_order=16), i, 17) for i in range(16)]),
     [0, 0, 1, 0, 0, 0, 0, 0], [-4, -4, 4, 4, 4, 4, -4, -4]),
    (17, 8, 5, 16, find_primitive_root(modulus=17, root_order=16),
     pow(find_primitive_root(modulus=17, root_order=16), 15, 17),
     [pow(find_primitive_root(modulus=17, root_order=16), i, 17) for i in range(16)],
     bit_reverse_copy([pow(find_primitive_root(modulus=17, root_order=16), i, 17) for i in range(16)]),
     [0, 0, 0, 1, 0, 0, 0, 0], [-2, 2, -2, 2, -8, 8, -8, 8]),
    (17, 8, 5, 16, find_primitive_root(modulus=17, root_order=16),
     pow(find_primitive_root(modulus=17, root_order=16), 15, 17),
     [pow(find_primitive_root(modulus=17, root_order=16), i, 17) for i in range(16)],
     bit_reverse_copy([pow(find_primitive_root(modulus=17, root_order=16), i, 17) for i in range(16)]),
     [0, 0, 0, 0, 1, 0, 0, 0], [-1, -1, -1, -1, 1, 1, 1, 1]),
    (17, 8, 5, 16, find_primitive_root(modulus=17, root_order=16),
     pow(find_primitive_root(modulus=17, root_order=16), 15, 17),
     [pow(find_primitive_root(modulus=17, root_order=16), i, 17) for i in range(16)],
     bit_reverse_copy([pow(find_primitive_root(modulus=17, root_order=16), i, 17) for i in range(16)]),
     [0, 0, 0, 0, 0, 1, 0, 0], [8, -8, -8, 8, -2, 2, 2, -2]),
    (17, 8, 5, 16, find_primitive_root(modulus=17, root_order=16),
     pow(find_primitive_root(modulus=17, root_order=16), 15, 17),
     [pow(find_primitive_root(modulus=17, root_order=16), i, 17) for i in range(16)],
     bit_reverse_copy([pow(find_primitive_root(modulus=17, root_order=16), i, 17) for i in range(16)]),
     [0, 0, 0, 0, 0, 0, 1, 0], [4, 4, -4, -4, 4, 4, -4, -4]),
    (17, 8, 5, 16, find_primitive_root(modulus=17, root_order=16),
     pow(find_primitive_root(modulus=17, root_order=16), 15, 17),
     [pow(find_primitive_root(modulus=17, root_order=16), i, 17) for i in range(16)],
     bit_reverse_copy([pow(find_primitive_root(modulus=17, root_order=16), i, 17) for i in range(16)]),
     [0, 0, 0, 0, 0, 0, 0, 1], [2, -2, 2, -2, -8, 8, -8, 8]),
]
GENTLEMAN_SANDE_INTT_TEST_DATA = [
    (17, 8, 5, 16, find_primitive_root(modulus=17, root_order=16),
     pow(find_primitive_root(modulus=17, root_order=16), 15, 17),
     [pow(pow(find_primitive_root(modulus=17, root_order=16), 15, 17), i, 17) for i in range(16)],
     bit_reverse_copy([pow(pow(find_primitive_root(modulus=17, root_order=16), 15, 17), i, 17) for i in range(16)]),
     [1, 0, 0, 0, 0, 0, 0, 0], [-2, -4, -8, 1, 2, 4, 8, -1]),
    (17, 8, 5, 16, find_primitive_root(modulus=17, root_order=16),
     pow(find_primitive_root(modulus=17, root_order=16), 15, 17),
     [pow(pow(find_primitive_root(modulus=17, root_order=16), 15, 17), i, 17) for i in range(16)],
     bit_reverse_copy([pow(pow(find_primitive_root(modulus=17, root_order=16), 15, 17), i, 17) for i in range(16)]),
     [0, 1, 0, 0, 0, 0, 0, 0], [-2, 4, -8, -1, 2, -4, 8, 1]),
    (17, 8, 5, 16, find_primitive_root(modulus=17, root_order=16),
     pow(find_primitive_root(modulus=17, root_order=16), 15, 17),
     [pow(pow(find_primitive_root(modulus=17, root_order=16), 15, 17), i, 17) for i in range(16)],
     bit_reverse_copy([pow(pow(find_primitive_root(modulus=17, root_order=16), 15, 17), i, 17) for i in range(16)]),
     [0, 0, 1, 0, 0, 0, 0, 0], [-2, 4, 8, 1, 2, -4, -8, -1]),
    (17, 8, 5, 16, find_primitive_root(modulus=17, root_order=16),
     pow(find_primitive_root(modulus=17, root_order=16), 15, 17),
     [pow(pow(find_primitive_root(modulus=17, root_order=16), 15, 17), i, 17) for i in range(16)],
     bit_reverse_copy([pow(pow(find_primitive_root(modulus=17, root_order=16), 15, 17), i, 17) for i in range(16)]),
     [0, 0, 0, 1, 0, 0, 0, 0], [-2, -4, 8, -1, 2, 4, -8, 1]),
    (17, 8, 5, 16, find_primitive_root(modulus=17, root_order=16),
     pow(find_primitive_root(modulus=17, root_order=16), 15, 17),
     [pow(pow(find_primitive_root(modulus=17, root_order=16), 15, 17), i, 17) for i in range(16)],
     bit_reverse_copy([pow(pow(find_primitive_root(modulus=17, root_order=16), 15, 17), i, 17) for i in range(16)]),
     [0, 0, 0, 0, 1, 0, 0, 0], [-2, 1, 8, -4, -2, 1, 8, -4]),
    (17, 8, 5, 16, find_primitive_root(modulus=17, root_order=16),
     pow(find_primitive_root(modulus=17, root_order=16), 15, 17),
     [pow(pow(find_primitive_root(modulus=17, root_order=16), 15, 17), i, 17) for i in range(16)],
     bit_reverse_copy([pow(pow(find_primitive_root(modulus=17, root_order=16), 15, 17), i, 17) for i in range(16)]),
     [0, 0, 0, 0, 0, 1, 0, 0], [-2, -1, 8, 4, -2, -1, 8, 4]),
    (17, 8, 5, 16, find_primitive_root(modulus=17, root_order=16),
     pow(find_primitive_root(modulus=17, root_order=16), 15, 17),
     [pow(pow(find_primitive_root(modulus=17, root_order=16), 15, 17), i, 17) for i in range(16)],
     bit_reverse_copy([pow(pow(find_primitive_root(modulus=17, root_order=16), 15, 17), i, 17) for i in range(16)]),
     [0, 0, 0, 0, 0, 0, 1, 0], [-2, -1, -8, -4, -2, -1, -8, -4]),
    (17, 8, 5, 16, find_primitive_root(modulus=17, root_order=16),
     pow(find_primitive_root(modulus=17, root_order=16), 15, 17),
     [pow(pow(find_primitive_root(modulus=17, root_order=16), 15, 17), i, 17) for i in range(16)],
     bit_reverse_copy([pow(pow(find_primitive_root(modulus=17, root_order=16), 15, 17), i, 17) for i in range(16)]),
     [0, 0, 0, 0, 0, 0, 0, 1], [-2, 1, -8, 4, -2, 1, -8, 4]),
]
NTT_POLY_MULT_TEST_DATA = [
    (17, 16, find_primitive_root(modulus=17, root_order=16),
     pow(find_primitive_root(modulus=17, root_order=16), 15, 17)),
]


@pytest.mark.parametrize('number, expected_result', IS_ODD_PRIME_TEST_DATA)
def test_is_odd_prime(number, expected_result):
    assert is_odd_prime(val=number) == expected_result


@pytest.mark.parametrize('number, expected_result', IS_POW_TWO_GEQ_TWO_TEST_DATA)
def test_is_pow_two_geq_two(number, expected_result):
    assert is_pow_two_geq_two(val=number) == expected_result


@pytest.mark.parametrize('modulus, root_order, expected_result', HAS_PRIMITIVE_ROOT_OF_UNITY_TEST_DATA)
def test_has_primitive_root_of_unity(modulus, root_order, expected_result):
    assert has_primitive_root_of_unity(modulus=modulus, root_order=root_order) == expected_result


@pytest.mark.parametrize('some_list, expected_result', BIT_REVERSE_COPY_TEST_DATA)
def test_bit_reverse_copy(some_list, expected_result):
    assert bit_reverse_copy(val=some_list) == expected_result
    assert bit_reverse_copy(val=expected_result) == some_list


@pytest.mark.parametrize('val, modulus, halfmod, logmod, expected_value', CENT_TEST_DATA)
def test_cent(val, modulus, halfmod, logmod, expected_value):
    assert halfmod == modulus//2
    assert 2**(logmod-1) < modulus <= 2**logmod
    assert cent(val=val, modulus=modulus, halfmod=halfmod, logmod=logmod) == expected_value


@pytest.mark.parametrize('val, modulus, root_order, expected_value', IS_ROOT_OF_UNITY_TEST_DATA)
def test_is_root_of_unity(val, modulus, root_order, expected_value):
    assert is_root_of_unity(purported_root=val, modulus=modulus, root_order=root_order) == expected_value


@pytest.mark.parametrize('val, modulus, root_order, expected_value', IS_PRIMITIVE_ROOT_OF_UNITY_TEST_DATA)
def test_is_primitive_root_of_unity(val, modulus, root_order, expected_value):
    assert is_primitive_root(purported_root=val, modulus=modulus, root_order=root_order) == expected_value


@pytest.mark.parametrize('modulus, root_order, expected_value', FIND_PRIMITIVE_ROOT_TEST_DATA)
def test_find_primitive_root(modulus, root_order, expected_value):
    assert find_primitive_root(modulus=modulus, root_order=root_order) == expected_value


@pytest.mark.parametrize(
    'modulus, halfmod, logmod, root_order, root, root_inv, root_powers, bit_rev_root_powers, val, expected_val',
    COOLEY_TUKEY_NTT_TEST_DATA)
def test_cooley_tukey(modulus, halfmod, logmod, root_order, root, root_inv, root_powers, bit_rev_root_powers, val,
                      expected_val):
    observed_val = cooley_tukey_ntt(val=val, modulus=modulus, root_order=root_order,
                                    bit_rev_root_powers=bit_rev_root_powers, halfmod=halfmod, logmod=logmod)
    assert all((x-y) % modulus == 0 for x, y in zip(observed_val, expected_val))


@pytest.mark.parametrize(
    'modulus, halfmod, logmod, root_order, root, root_inv, root_inv_powers, bit_rev_root_inv_powers, val, expected_val',
    GENTLEMAN_SANDE_INTT_TEST_DATA)
def test_gentleman_sande(modulus, halfmod, logmod, root_order, root, root_inv, root_inv_powers, bit_rev_root_inv_powers,
                         val, expected_val):
    observed_val = gentleman_sande_intt(val=val, modulus=modulus, root_order=root_order,
                                        bit_rev_inv_root_powers=bit_rev_root_inv_powers, halfmod=halfmod, logmod=logmod)
    assert all((x-y) % modulus == 0 for x, y in zip(observed_val, expected_val))


@pytest.fixture(params=NTT_POLY_MULT_TEST_DATA)
def some_ran_polys_and_prod(request):
    modulus, root_order, root, root_inv = request.param
    degree = root_order//2
    left_factors: List[List[int]] = []
    right_factors: List[List[int]] = []
    expected_product_by_foil: List[List[int]] = []
    while len(left_factors) < SAMPLE_SIZE:
        new_left_factor = [randbelow(modulus) for _ in range(degree)]
        left_factors += [new_left_factor]

        new_right_factor = [randbelow(modulus) for _ in range(degree)]
        right_factors += [new_right_factor]

        new_expected_product_by_foil = [0 for _ in range(root_order)]
        for i, a in enumerate(new_left_factor):
            for j, b in enumerate(new_right_factor):
                new_expected_product_by_foil[i+j] += a*b % modulus
        new_expected_product_by_foil = [(x - y) % modulus for x, y in zip(new_expected_product_by_foil[:degree],
                                                                          new_expected_product_by_foil[degree:])]
        expected_product_by_foil += [new_expected_product_by_foil]

        assert len(expected_product_by_foil[-1]) == degree
    return modulus, root_order, root, root_inv, left_factors, right_factors, expected_product_by_foil


def test_ntt_poly_mult(some_ran_polys_and_prod):
    modulus, root_order, root, inv_root, left_factors, right_factors, expected_vals = some_ran_polys_and_prod
    for left_factor, right_factor, expected_val in zip(left_factors, right_factors, expected_vals):
        halfmod = modulus//2
        logmod = modulus.bit_length()
        degree = len(left_factor)
        root_powers = [pow(base=root, exp=i, mod=modulus) for i in range(root_order)]
        inv_root_powers = [pow(base=inv_root, exp=i, mod=modulus) for i in range(root_order)]
        brv_root_powers = bit_reverse_copy(root_powers)
        brv_inv_root_powers = bit_reverse_copy(inv_root_powers)
        observed_val = ntt_poly_mult(f=left_factor, g=right_factor, modulus=modulus, halfmod=halfmod, logmod=logmod,
                                     degree=degree, root_order=2*degree, root=root, inv_root=inv_root,
                                     brv_root_powers=brv_root_powers)
        assert all((x-y) % modulus == 0 for x, y in zip(observed_val, expected_val))
