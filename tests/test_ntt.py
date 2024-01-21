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
from copy import deepcopy

SAMPLE_SIZE: int = 2**2
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
Q_MAX: int = 2**23
PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE: List[Tuple[int, int]] = []
for log2d in range(LOG2_D_MIN, LOG2_D_MAX + 1):
    tmp_d: int = 1 << log2d
    tmp_q: int = 2 * tmp_d + 1
    while tmp_q < Q_MAX:
        while not is_odd_prime(tmp_q) and tmp_q < Q_MAX:
            tmp_q += 2 * tmp_d
        if is_odd_prime(tmp_q) and tmp_q < Q_MAX:
            PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE.append((tmp_d, tmp_q))
            find_primitive_root(modulus=tmp_q, root_order=2 * tmp_d)  # precompute and cache
            tmp_q *= 8
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
TEST_CASES_Q17D8 = [
    (17, 8, 17//2, (17).bit_length())  # 0 modulus, 1 degree, 2 halfmod, 3 logmod
]
TEST_CASES_Q17D8 = [
    t + tuple([find_primitive_root(modulus=t[0], root_order=2*t[1])])  # 4 root
    for t in TEST_CASES_Q17D8
]
TEST_CASES_Q17D8 = [
    t + tuple([bit_reverse_copy([pow(base=t[-1], exp=i, mod=t[0]) for i in range(t[1])])])  # 5 brv_root_powers
    for t in TEST_CASES_Q17D8
]
TEST_CASES_Q17D8 = [
    t + tuple([pow(base=t[4], exp=t[0]-2, mod=t[0])])  # 6 inv_root
    for t in TEST_CASES_Q17D8
]
TEST_CASES_Q17D8 = [
    t + tuple([bit_reverse_copy([pow(base=t[-1], exp=i, mod=t[0]) for i in range(t[1])])])  # 7 brv_inv_root_powers
    for t in TEST_CASES_Q17D8
]
TEST_RAN_MUL = [
    t + tuple([
        [randbelow(t[0]) for _ in range(t[1])],
        [randbelow(t[0]) for _ in range(t[1])]
    ])
    for t in TEST_CASES_Q17D8
    for _ in range(SAMPLE_SIZE)
]
TEST_RAN_MUL_CASES = []
for t in TEST_RAN_MUL:
    f = t[-2]
    g = t[-1]
    fg = [0 for _ in range(2*t[1])]
    for i, a in enumerate(f):
        for j, b in enumerate(g):
            fg[i+j] += a*b
    fg = [cent(val=x-y, modulus=t[0], halfmod=t[2], logmod=t[3]) for x, y in zip(fg[:t[1]], fg[t[1]:])]
    TEST_RAN_MUL_CASES += [t + tuple([
        fg,
        cooley_tukey_ntt(val=deepcopy(f), modulus=t[0], root_order=2*t[1], bit_rev_root_powers=t[5], halfmod=t[2], logmod=t[3]),
        cooley_tukey_ntt(val=deepcopy(g), modulus=t[0], root_order=2*t[1], bit_rev_root_powers=t[5], halfmod=t[2], logmod=t[3]),
        cooley_tukey_ntt(val=deepcopy(fg), modulus=t[0], root_order=2*t[1], bit_rev_root_powers=t[5], halfmod=t[2], logmod=t[3]),
    ])]
TEST_RAN_MUL_CASES = [t + tuple([
    [cent(val=x*y, modulus=t[0], halfmod=t[2], logmod=t[3]) for x, y in zip(t[-3], t[-2])]
]) for t in TEST_RAN_MUL_CASES]
TEST_RAN_MUL_CASES = [t + tuple([
    gentleman_sande_intt(val=deepcopy(t[-1]), modulus=t[0], root_order=2*t[1], bit_rev_inv_root_powers=t[7], halfmod=t[2], logmod=t[3]),
]) for t in TEST_RAN_MUL_CASES]
TEST_RAN_MUL_CASES = [t + tuple([
    [cent(val=x-y, modulus=t[0], halfmod=t[2], logmod=t[3]) for x, y in zip(t[-1], t[-6])]
]) for t in TEST_RAN_MUL_CASES]


@pytest.mark.parametrize('modulus,degree,halfmod,logmod,root,brv_root_powers,inv_root,brv_inv_root_powers,' +
                         'f,g,hand_fg,ntt_f,ntt_g,ntt_hand_fg,hadamard_fhat_ghat,intt_hadamard,diff', TEST_RAN_MUL_CASES)
def test_ran_mul(modulus,degree,halfmod,logmod,root,brv_root_powers,inv_root,brv_inv_root_powers,f,g,hand_fg,ntt_f,ntt_g,ntt_hand_fg,hadamard_fhat_ghat,intt_hadamard,diff):
    print()
    assert all(x % modulus == 0 for x in diff)
    assert all((tomato - carrot) % modulus == 0 for tomato, carrot in zip(hadamard_fhat_ghat,ntt_hand_fg))  # NTT(f) star NTT(g) = NTT(f * g)
    assert all((tomato - carrot) % modulus == 0 for tomato, carrot in zip(intt_hadamard, hand_fg))  # INTT(NTT(f) star NTT(g)) = f * g



TEST_MONOMIALS = [
    [1, 0, 0, 0, 0, 0, 0, 0,], [0, 1, 0, 0, 0, 0, 0, 0,], [0, 0, 1, 0, 0, 0, 0, 0,], [0, 0, 0, 1, 0, 0, 0, 0,],
    [0, 0, 0, 0, 1, 0, 0, 0,], [0, 0, 0, 0, 0, 1, 0, 0,], [0, 0, 0, 0, 0, 0, 1, 0,], [0, 0, 0, 0, 0, 0, 0, 1,]]
TEST_FORWARD_MONOMIALS = [
    [1, 1, 1, 1, 1, 1, 1, 1], [3, -3, 5, -5, -7, 7, -6, 6], [-8, -8, 8, 8, -2, -2, 2, 2], [-7, 7, 6, -6, -3, 3, 5, -5],
    [-4, -4, -4, -4, 4, 4, 4, 4], [5, -5, -3, 3, 6, -6, -7, 7], [-2, -2, 2, 2, -8, -8, 8, 8], [-6, 6, -7, 7, 5, -5, 3, -3]]  # NTT(monomials)
TEST_BACKWARD_MONOMIALS = [
    [-2, 5, -4, -7, -8, 3, 1, 6], [-2, -5, -4, 7, -8, -3, 1, -6], [-2, 3, 4, -6, -8, -5, -1, -7], [-2, -3, 4, 6, -8, 5, -1, 7],
    [-2, -7, 1, -5, 8, -6, -4, 3], [-2, 7, 1, 5, 8, 6, -4, -3], [-2, 6, -1, 3, 8, -7, 4, 5], [-2, -6, -1, -3, 8, 7, 4, -5]]  # INTT(monomials)
TEST_CASES_Q17D8 = [
    t + tuple([
        TEST_MONOMIALS[i],  # a unit vector
        TEST_FORWARD_MONOMIALS[i],  # the hand-computed forward transform of a
        TEST_BACKWARD_MONOMIALS[i], # the hand-computed inverse transform of a
    ])
    for t in TEST_CASES_Q17D8
    for i in range(len(TEST_MONOMIALS))
]
TEST_CASES_Q17D8 = [
    t + tuple([
        cooley_tukey_ntt(val=deepcopy(t[-3]), modulus=t[0], root_order=2*t[1], bit_rev_root_powers=t[5], halfmod=t[2], logmod=t[3]),  # forward transform of a
        gentleman_sande_intt(val=deepcopy(t[-3]), modulus=t[0], root_order=2*t[1], bit_rev_inv_root_powers=t[7], halfmod=t[2], logmod=t[3]),  # inverse transform of a
        gentleman_sande_intt(val=deepcopy(t[-2]), modulus=t[0], root_order=2*t[1], bit_rev_inv_root_powers=t[7], halfmod=t[2], logmod=t[3]),  # inverse transform of hand-computed forward transform
        cooley_tukey_ntt(val=deepcopy(t[-1]), modulus=t[0], root_order=2*t[1], bit_rev_root_powers=t[5], halfmod=t[2], logmod=t[3]),   # forward transform of hand-computed inverse transform
    ])
    for t in TEST_CASES_Q17D8
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


@pytest.mark.parametrize('modulus,degree,halfmod,logmod,root,brv_root_powers,inv_root,brv_inv_root_powers,' +
                         'input_vector,hand_ntt_input_vector,hand_intt_input_vector,ntt_input_vector,' +
                         'intt_input_vector,intt_hand_ntt_input_vector,ntt_hand_intt_input_vector', TEST_CASES_Q17D8)
def test_q17d18(modulus, degree, halfmod, logmod, root, brv_root_powers, inv_root, brv_inv_root_powers,
                input_vector, hand_ntt_input_vector, hand_intt_input_vector, ntt_input_vector,
                intt_input_vector, intt_hand_ntt_input_vector, ntt_hand_intt_input_vector):
    assert all((tomato - carrot) % modulus == 0 for tomato, carrot in zip(hand_ntt_input_vector, ntt_input_vector))
    assert all((tomato - carrot) % modulus == 0 for tomato, carrot in zip(hand_intt_input_vector, intt_input_vector))
    assert all((tomato - carrot) % modulus == 0 for tomato, carrot in zip(input_vector, intt_hand_ntt_input_vector))
    assert all((tomato - carrot) % modulus == 0 for tomato, carrot in zip(input_vector, ntt_hand_intt_input_vector))


@pytest.mark.parametrize('modulus,degree,halfmod,logmod,root,brv_root_powers,inv_root,brv_inv_root_powers,' +
                         'f,g,hand_fg,ntt_f,ntt_g,ntt_hand_fg,hadamard_fhat_ghat,intt_hadamard,diff', TEST_RAN_MUL_CASES)
def test_ntt_poly_mult(modulus,degree,halfmod,logmod,root,brv_root_powers,inv_root,brv_inv_root_powers,f,g,hand_fg,ntt_f,ntt_g,ntt_hand_fg,hadamard_fhat_ghat,intt_hadamard,diff):
    observed_val = ntt_poly_mult(f=f, g=g, modulus=modulus, halfmod=halfmod, logmod=logmod,
                                 degree=degree, root_order=2*degree, root=root, inv_root=inv_root,
                                 brv_root_powers=brv_root_powers)
    assert all((a-b) % modulus == 0 for a, b in zip(observed_val, hand_fg))

