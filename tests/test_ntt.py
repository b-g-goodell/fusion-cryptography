import pytest
from typing import List, Tuple
from secrets import randbelow
from algebra.ntt import (
    _is_odd_prime,
    _is_pos_pow_two,
    is_ntt_friendly,
    _brv_indices,
    brv_copy,
    _brv_root_and_inv_root_powers,
    cent,
    _is_rou,
    _is_prou,
    find_prou,
    _cooley_tukey_ntt,
    _gentleman_sande_intt,
    _ntt_poly_mult,
    ntt
)
from copy import deepcopy

SAMPLE_SIZE: int = 2 ** 10
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
IS_POS_POW_TWO_TEST_CASES = [
    (3, False),
    (2, True),
    (1, True),
    (4, True),
    (8, True),
    (7, False)
]
LOG2_D_MIN: int = 2
LOG2_D_MAX: int = 6
Q_MAX: int = 2 ** 23
PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE: List[Tuple[int, int]] = []
for log2d in range(LOG2_D_MIN, LOG2_D_MAX + 1):
    tmp_d: int = 1 << log2d
    tmp_q: int = 2 * tmp_d + 1
    while tmp_q < Q_MAX:
        while not _is_odd_prime(tmp_q) and tmp_q < Q_MAX:
            tmp_q += 2 * tmp_d
        if _is_odd_prime(tmp_q) and tmp_q < Q_MAX:
            PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE.append((tmp_d, tmp_q))
            find_prou(mod=tmp_q, deg=tmp_d)
            tmp_q *= 2
            tmp_q -= (tmp_q - 1) % (2 * tmp_d)
            assert (tmp_q - 1) % (2 * tmp_d) == 0
HAS_PRIMITIVE_ROOT_OF_UNITY_TEST_DATA = [
    (17, 2, True),
    (17, 4, True),
    (17, 8, True),
    (17, 16, False),
    (17, 3, False),
    (17, 5, False),
    (17, 6, False),
    (17, 7, False),
    (257, 128, True),
    (257, 256, False)
]
BIT_REVERSE_COPY_TEST_DATA = [
    (4, [0, 2, 1, 3], _brv_indices(val=4)),
    (8, [0, 4, 2, 6, 1, 5, 3, 7], _brv_indices(val=8))
]
CENT_TEST_DATA = [
    (0, 17, 0),
    (1, 17, 1),
    (2, 17, 2),
    (3, 17, 3),
    (4, 17, 4),
    (5, 17, 5),
    (6, 17, 6),
    (7, 17, 7),
    (8, 17, 8),
    (9, 17, -8),
    (10, 17, -7),
    (11, 17, -6),
    (12, 17, -5),
    (13, 17, -4),
    (14, 17, -3),
    (15, 17, -2),
    (16, 17, -1),
]
IS_ROOT_OF_UNITY_TEST_DATA = [
    (3, 17, 8, True),
    (2, 17, 8, True),
]
IS_PRIMITIVE_ROOT_OF_UNITY_TEST_DATA = [
    (3, 17, 8, True),
    (2, 17, 8, False),
]
FIND_PRIMITIVE_ROOT_TEST_DATA = [
    (17, 8, 3),
    (257, 128, 3)
]
COOLEY_TUKEY_NTT_TEST_DATA = [
    (17, 8)
]
COOLEY_TUKEY_NTT_TEST_DATA = [
    t + tuple([
        find_prou(mod=t[0], deg=t[1])
    ])
    for t in COOLEY_TUKEY_NTT_TEST_DATA
]
COOLEY_TUKEY_NTT_TEST_DATA = [
    t + tuple([
        pow(base=t[-1], exp=t[0]-2, mod=t[0]),
        brv_copy(val=_brv_root_and_inv_root_powers(deg=t[1], mod=t[0])[0]),
        _brv_root_and_inv_root_powers(deg=t[1], mod=t[0])[0],
        brv_copy(val=_brv_root_and_inv_root_powers(deg=t[1], mod=t[0])[1]),
        _brv_root_and_inv_root_powers(deg=t[1], mod=t[0])[1]
    ])
    for t in COOLEY_TUKEY_NTT_TEST_DATA
]
assert len(COOLEY_TUKEY_NTT_TEST_DATA) == 1
COOLEY_TUKEY_NTT_TEST_DATA = [
    COOLEY_TUKEY_NTT_TEST_DATA[0] + tuple([[1, 0, 0, 0, 0, 0, 0, 0], [1, 1, 1, 1, 1, 1, 1, 1]]),
    COOLEY_TUKEY_NTT_TEST_DATA[0] + tuple([[0, 1, 0, 0, 0, 0, 0, 0], [3, -3, 5, -5, -7, 7, -6, 6]]),
    COOLEY_TUKEY_NTT_TEST_DATA[0] + tuple([[0, 0, 1, 0, 0, 0, 0, 0], [-8, -8, 8, 8, -2, -2, 2, 2]]),
    COOLEY_TUKEY_NTT_TEST_DATA[0] + tuple([[0, 0, 0, 1, 0, 0, 0, 0], [-7, 7, 6, -6, -3, 3, 5, -5]]),
    COOLEY_TUKEY_NTT_TEST_DATA[0] + tuple([[0, 0, 0, 0, 1, 0, 0, 0], [-4, -4, -4, -4, 4, 4, 4, 4]]),
    COOLEY_TUKEY_NTT_TEST_DATA[0] + tuple([[0, 0, 0, 0, 0, 1, 0, 0], [5, -5, -3, 3, 6, -6, -7, 7]]),
    COOLEY_TUKEY_NTT_TEST_DATA[0] + tuple([[0, 0, 0, 0, 0, 0, 1, 0], [-2, -2, 2, 2, -8, -8, 8, 8]]),
    COOLEY_TUKEY_NTT_TEST_DATA[0] + tuple([[0, 0, 0, 0, 0, 0, 0, 1], [-6, 6, -7, 7, 5, -5, 3, -3]]),
]
TEST_NTT_DATA = [
    (17, [1, 0, 0, 0, 0, 0, 0, 0], False, [1, 1, 1, 1, 1, 1, 1, 1],
     ntt(val=[1, 0, 0, 0, 0, 0, 0, 0], mod=17, inv_flag=False)),
    (17, [0, 1, 0, 0, 0, 0, 0, 0], False, [3,-3,5,-5,-7,7,-6,6],
     ntt(val=[0, 1, 0, 0, 0, 0, 0, 0], mod=17, inv_flag=False)),
    (17, [0, 0, 1, 0, 0, 0, 0, 0], False, [-8,-8,8,8,-2,-2,2,2],
     ntt(val=[0, 0, 1, 0, 0, 0, 0, 0], mod=17, inv_flag=False)),
    (17, [0, 0, 0, 1, 0, 0, 0, 0], False, [-7,7,6,-6,-3,3,5,-5],
     ntt(val=[0, 0, 0, 1, 0, 0, 0, 0], mod=17, inv_flag=False)),
    (17, [0, 0, 0, 0, 1, 0, 0, 0], False, [-4,-4,-4,-4,4,4,4,4],
     ntt(val=[0, 0, 0, 0, 1, 0, 0, 0], mod=17, inv_flag=False)),
    (17, [0, 0, 0, 0, 0, 1, 0, 0], False, [5, -5, -3, 3, 6, -6, -7, 7],
     ntt(val=[0, 0, 0, 0, 0, 1, 0, 0], mod=17, inv_flag=False)),
    (17, [0, 0, 0, 0, 0, 0, 1, 0], False, [-2, -2, 2, 2, -8, -8, 8, 8],
     ntt(val=[0, 0, 0, 0, 0, 0, 1, 0], mod=17, inv_flag=False)),
    (17, [0, 0, 0, 0, 0, 0, 0, 1], False, [-6, 6, -7, 7, 5, -5, 3, -3],
     ntt(val=[0, 0, 0, 0, 0, 0, 0, 1], mod=17, inv_flag=False)),
    (17, [1, 1, 1, 1, 1, 1, 1, 1], True, [1, 0, 0, 0, 0, 0, 0, 0],
     ntt(val=[1, 1, 1, 1, 1, 1, 1, 1], mod=17, inv_flag=True)),
    (17, [3,-3,5,-5,-7,7,-6,6], True, [0, 1, 0, 0, 0, 0, 0, 0],
     ntt(val=[3, -3, 5, -5, -7, 7, -6, 6], mod=17, inv_flag=True)),
    (17, [-8,-8,8,8,-2,-2,2,2], True, [0, 0, 1, 0, 0, 0, 0, 0],
     ntt(val=[-8, -8, 8, 8, -2, -2, 2, 2], mod=17, inv_flag=True)),
    (17, [-7,7,6,-6,-3,3,5,-5], True, [0, 0, 0, 1, 0, 0, 0, 0],
     ntt(val=[-7, 7, 6, -6, -3, 3, 5, -5], mod=17, inv_flag=True)),
    (17, [-4,-4,-4,-4,4,4,4,4], True, [0, 0, 0, 0, 1, 0, 0, 0],
     ntt(val=[-4, -4, -4, -4, 4, 4, 4, 4], mod=17, inv_flag=True)),
    (17, [5, -5, -3, 3, 6, -6, -7, 7], True, [0, 0, 0, 0, 0, 1, 0, 0],
     ntt(val=[5, -5, -3, 3, 6, -6, -7, 7], mod=17, inv_flag=True)),
    (17, [-2, -2, 2, 2, -8, -8, 8, 8], True, [0, 0, 0, 0, 0, 0, 1, 0],
     ntt(val=[-2, -2, 2, 2, -8, -8, 8, 8], mod=17, inv_flag=True)),
    (17, [-6, 6, -7, 7, 5, -5, 3, -3], True, [0, 0, 0, 0, 0, 0, 0, 1],
     ntt(val=[-6, 6, -7, 7, 5, -5, 3, -3], mod=17, inv_flag=True)),
]

@pytest.mark.parametrize('number, expected_result', IS_ODD_PRIME_TEST_DATA)
def test_is_odd_prime(number, expected_result):
    assert _is_odd_prime(val=number) == expected_result


@pytest.mark.parametrize('modulus, degree, expected_result', HAS_PRIMITIVE_ROOT_OF_UNITY_TEST_DATA)
def test_has_prou(modulus, degree, expected_result):
    assert is_ntt_friendly(mod=modulus, deg=degree) == expected_result


@pytest.mark.parametrize('number, expected_result', IS_POS_POW_TWO_TEST_CASES)
def test_is_pos_pow_two(number, expected_result):
    assert _is_pos_pow_two(val=number) == expected_result


@pytest.mark.parametrize('list_len, expected_result, observed_result', BIT_REVERSE_COPY_TEST_DATA)
def test_bit_reverse_len_to_idxs(list_len, expected_result, observed_result):
    assert expected_result == observed_result


def test_bit_reverse_copy():
    pass


@pytest.mark.parametrize('val, modulus, degree, expected_value', IS_ROOT_OF_UNITY_TEST_DATA)
def test_is_root_of_unity(val, modulus, degree, expected_value):
    assert _is_rou(root=val, mod=modulus, deg=degree) == expected_value


@pytest.mark.parametrize('val, modulus, degree, expected_value', IS_PRIMITIVE_ROOT_OF_UNITY_TEST_DATA)
def test_is_primitive_root_of_unity(val, modulus, degree, expected_value):
    assert _is_prou(root=val, mod=modulus, deg=degree) == expected_value


@pytest.mark.parametrize('modulus, deg, expected_value', FIND_PRIMITIVE_ROOT_TEST_DATA)
def test_find_primitive_root(modulus, deg, expected_value):
    assert find_prou(mod=modulus, deg=deg) == expected_value


def test_is_root_inverse():
    pass


def test_is_ntt_valid():
    pass


@pytest.mark.parametrize('val, modulus, expected_value', CENT_TEST_DATA)
def test_cent(val, modulus, expected_value):
    assert cent(val=val, mod=modulus) == expected_value


@pytest.mark.parametrize(
    'modulus, deg, root, root_inv, root_powers, bit_rev_root_powers, inv_root_powers, bit_rev_inv_root_powers, val, expected_val',
    COOLEY_TUKEY_NTT_TEST_DATA)
def test_cooley_tukey(modulus, deg, root, root_inv, root_powers, bit_rev_root_powers, inv_root_powers, bit_rev_inv_root_powers, val, expected_val):
    dcval = deepcopy(val)
    observed_val = _cooley_tukey_ntt(val=dcval, mod=modulus, brv_powers=bit_rev_root_powers)
    assert all((x - y) % modulus == 0 for x, y in zip(observed_val, expected_val))


@pytest.mark.parametrize('deg, mod', PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE)
def test_random_poly_mul(deg, mod):
    assert is_ntt_friendly(mod=mod, deg=deg)
    root = find_prou(mod=mod, deg=deg)
    inv_root = pow(base=root, exp=mod - 2, mod=mod)
    root_powers = [pow(base=root, exp=i, mod=mod) for i in range(2 * deg + 1)]
    inv_root_powers = [pow(base=inv_root, exp=i, mod=mod) for i in range(2 * deg + 1)]
    brv_powers = brv_copy(root_powers[:deg])
    brv_inv_root_powers = brv_copy(inv_root_powers[:deg])

    for _ in range(SAMPLE_SIZE):
        left_factor = [randbelow(mod) for _ in range(deg)]
        right_factor = [randbelow(mod) for _ in range(deg)]
        exp_product_from_foil = [0 for _ in range(2 * deg)]
        for i, a in enumerate(left_factor):
            for j, b in enumerate(right_factor):
                exp_product_from_foil[i + j] += a * b % mod
        exp_product_from_foil = [(x - y) % mod for x, y in
                                 zip(exp_product_from_foil[:deg], exp_product_from_foil[deg:])]
        obs_prod_from_ntt_poly_mult = _ntt_poly_mult(f=left_factor, g=right_factor, mod=mod, brv_powers=brv_powers, brv_inv_root_powers=brv_inv_root_powers)
        assert all((x - y) % mod == 0 for x, y in zip(obs_prod_from_ntt_poly_mult, exp_product_from_foil))


# Integration test
@pytest.mark.parametrize(
    'modulus, deg, root, root_inv, root_powers, bit_rev_root_powers, inv_root_powers, bit_rev_inv_root_powers, val, expected_val',
    COOLEY_TUKEY_NTT_TEST_DATA)
def test_ct_gs_inverse(modulus, deg, root, root_inv, root_powers, bit_rev_root_powers, inv_root_powers, bit_rev_inv_root_powers, val, expected_val):
    x = deepcopy(val)
    y = deepcopy(x)
    _cooley_tukey_ntt(val=x, mod=modulus, brv_powers=bit_rev_root_powers)
    assert x != y
    inv_root_powers = [pow(base=root_inv, exp=i, mod=modulus) for i in range(deg)]
    brv_inv_root_powers = brv_copy(inv_root_powers)
    _gentleman_sande_intt(val=x, mod=modulus, brv_powers=brv_inv_root_powers)
    assert x == y


@pytest.mark.parametrize("modulus, val, inv_flag, expected_result, observed_result", TEST_NTT_DATA)
def test_ntt(modulus, val, inv_flag, expected_result, observed_result):
    assert observed_result == expected_result
    z = ntt(val=val, mod=modulus, inv_flag=inv_flag)
    assert all((x-y) % modulus == 0 for x, y in zip(z, observed_result))