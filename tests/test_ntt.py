import pytest
from typing import List, Tuple
from secrets import randbelow
from algebra.ntt import (
    _is_odd_prime,
    _is_pow_two_geq_two,
    _has_primitive_root_of_unity,
    _bit_reverse_copy,
    _check_modulus_halfmod_logmod,
    _cent,
    _is_root_of_unity,
    _is_primitive_root,
    _find_primitive_root,
    _check_ntt_and_intt,
    _cooley_tukey_ntt,
    _gentleman_sande_intt,
    _ntt_poly_mult,
    _ERR_MUST_BE_INT,
    _ERR_MUST_BE_LIST,
    _ERR_MUST_BE_BOOL,
    _ERR_MUST_BE_POW_TWO,
    _ERR_MUST_BE_GEQ_THREE,
    _ERR_MUST_BE_STRICTLY_POSITIVE,
    _ERR_MUST_BE_FLOOR_HALF,
    _ERR_MUST_BE_BIT_LENGTH,
    _ERR_MUST_HAVE_PRIM_ROU,
    _ERR_MUST_BE_ODD_PRIME,
    _ERR_MUST_BE_CORRECT_ROOT,
    _ERR_MUST_BE_SAME_LENGTH,
)
from copy import deepcopy

SAMPLE_SIZE: int = 2**2
IS_ODD_PRIME_TEST_DATA: List[Tuple[int, bool]] = [
    (3, True),
    (4, False),
    (5, True),
    (6, False),
    (29, True),
    (30, False),
    (1, False),
    (2, False),
]
LOG2_D_MIN: int = 2
LOG2_D_MAX: int = 6
Q_MAX: int = 2**23
PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE: List[Tuple[int, int]] = []
for log2d in range(LOG2_D_MIN, LOG2_D_MAX + 1):
    tmp_d: int = 1 << log2d
    tmp_q: int = 2 * tmp_d + 1
    while tmp_q < Q_MAX:
        while not _is_odd_prime(tmp_q) and tmp_q < Q_MAX:
            tmp_q += 2 * tmp_d
        if _is_odd_prime(tmp_q) and tmp_q < Q_MAX:
            PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE.append((tmp_d, tmp_q))
            _find_primitive_root(modulus=tmp_q, root_order=2 * tmp_d)  # precompute and cache
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
IS_POW_TWO_GEQ_TWO_TEST_DATA = [
    (3, False),
    (2, True),
    (1, False),
    (4, True),
    (8, True),
    (7, False)
]
BIT_REVERSE_COPY_TEST_DATA = [
    ([0, 1, 2, 3], [0, 2, 1, 3]),
    ([0, 1, 2, 3, 4, 5, 6, 7], [0, 4, 2, 6, 1, 5, 3, 7])
]
CHECK_MODULUS_HALFMOD_LOGMOD_TEST_DATA = [
    (x, x // 2, x.bit_length())  # adjusted logmod to be x.bit_length() - 1
    for x in range(3, 100)
]
CHECK_MODULUS_HALFMOD_LOGMOD_ERRORS = [
    ("not an int", 3, 3, TypeError, _ERR_MUST_BE_INT),
    (6, "not an int", 3, TypeError, _ERR_MUST_BE_INT),
    (6, 3, "not an int", TypeError, _ERR_MUST_BE_INT),
    (2, 1, 2, ValueError, _ERR_MUST_BE_GEQ_THREE),
    (3, 0, 2, ValueError, _ERR_MUST_BE_STRICTLY_POSITIVE),
    (3, 2, 2, ValueError, _ERR_MUST_BE_FLOOR_HALF),
    (6, 2, 0, ValueError, _ERR_MUST_BE_STRICTLY_POSITIVE),
    (6, 3, 2, ValueError, _ERR_MUST_BE_BIT_LENGTH)  # adjust the error message to match the actual output
    # The above line has been wrapped for readability, but should be a single line in the code.
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
CHECK_NTT_AND_INTT_VALID_DATA = [
    (list(range(8)), 17, 16)   #  val, modulus, root_order = 2*degree
]
CHECK_NTT_AND_INTT_VALID_DATA = [
    t + tuple([
        False,
        _find_primitive_root(modulus=t[1], root_order=t[2]),  # root
    ])
    for t in CHECK_NTT_AND_INTT_VALID_DATA
] + [
    t + tuple([
        True,
        pow(base=_find_primitive_root(modulus=t[1], root_order=t[2]), exp=t[1] - 2, mod=t[1]),
    ])
    for t in CHECK_NTT_AND_INTT_VALID_DATA
]
CHECK_NTT_AND_INTT_VALID_DATA = [
    t + tuple([
        _bit_reverse_copy([pow(base=t[-1], exp=i, mod=t[1]) for i in range(t[2] // 2)]),  # brv_powers
        t[1]//2,  # halfmod
        t[1].bit_length()  # logmod
    ])
    for t in CHECK_NTT_AND_INTT_VALID_DATA
]
CHECK_NTT_AND_INTT_ERRORS = [
    ("not a list of ints", 17, 16, False, _bit_reverse_copy([pow(base=_find_primitive_root(modulus=17, root_order=16), exp=i, mod=17) for i in range(8)]), 8, 5, TypeError, _ERR_MUST_BE_LIST),
    (list("not a list of ints"), 17, 16, False, _bit_reverse_copy([pow(base=_find_primitive_root(modulus=17, root_order=16), exp=i, mod=17) for i in range(8)]), 8, 5, TypeError, _ERR_MUST_BE_INT),
    (list(range(8)), "not an int", 16, False, _bit_reverse_copy([pow(base=_find_primitive_root(modulus=17, root_order=16), exp=i, mod=17) for i in range(8)]), 8, 5, TypeError, _ERR_MUST_BE_INT),
    (list(range(8)), 17, "not an int", False, _bit_reverse_copy([pow(base=_find_primitive_root(modulus=17, root_order=16), exp=i, mod=17) for i in range(8)]), 8, 5, TypeError, _ERR_MUST_BE_INT),
    (list(range(8)), 17, 16, "not a bool", _bit_reverse_copy([pow(base=_find_primitive_root(modulus=17, root_order=16), exp=i, mod=17) for i in range(8)]), 8, 5, TypeError, _ERR_MUST_BE_BOOL),
    (list(range(8)), 17, 16, False, "not a list of ints", 8, 5, TypeError, _ERR_MUST_BE_LIST),
    (list(range(8)), 17, 16, False, list("not a list of ints"), 8, 5, TypeError, _ERR_MUST_BE_INT),
    (list(range(8)), 17, 16, False, _bit_reverse_copy([pow(base=_find_primitive_root(modulus=17, root_order=16), exp=i, mod=17) for i in range(8)]), "not an int", 5, TypeError, _ERR_MUST_BE_INT),
    (list(range(8)), 17, 16, False, _bit_reverse_copy([pow(base=_find_primitive_root(modulus=17, root_order=16), exp=i, mod=17) for i in range(8)]), 8, "not an int", TypeError, _ERR_MUST_BE_INT),
    (list(range(8)), 18, 16, False, _bit_reverse_copy([pow(base=_find_primitive_root(modulus=17, root_order=16), exp=i, mod=17) for i in range(8)]), 8, 5, ValueError, _ERR_MUST_BE_ODD_PRIME),
    (list(range(16)), 17, 32, False, _bit_reverse_copy([pow(base=_find_primitive_root(modulus=17, root_order=16), exp=i, mod=17) for i in range(16)]), 8, 5, ValueError, _ERR_MUST_HAVE_PRIM_ROU),
    (list(range(15)), 17, 16, False, _bit_reverse_copy([pow(base=_find_primitive_root(modulus=17, root_order=16), exp=i, mod=17) for i in range(8)]), 8, 5, ValueError, _ERR_MUST_BE_POW_TWO),
    (list(range(8)), 17, 8, False, _bit_reverse_copy([pow(base=_find_primitive_root(modulus=17, root_order=8), exp=i, mod=17) for i in range(8)]), 8, 5, ValueError, _ERR_MUST_BE_FLOOR_HALF),
    (list(range(8)), 17, 16, False, _bit_reverse_copy([pow(base=_find_primitive_root(modulus=17, root_order=16), exp=i, mod=17) for i in range(8)]), 7, 5, ValueError, _ERR_MUST_BE_FLOOR_HALF),
    (list(range(8)), 17, 16, False, _bit_reverse_copy([pow(base=_find_primitive_root(modulus=17, root_order=16), exp=i, mod=17) for i in range(8)]), 8, 4, ValueError, _ERR_MUST_BE_BIT_LENGTH),
    (list(range(8)), 17, 16, False, _bit_reverse_copy([1 + pow(base=_find_primitive_root(modulus=17, root_order=16), exp=i, mod=17) for i in range(8)]), 8, 5, ValueError, _ERR_MUST_BE_CORRECT_ROOT),
    (list(range(8)), 17, 16, False, _bit_reverse_copy([pow(base=pow(base=_find_primitive_root(modulus=17, root_order=16), exp=15, mod=17), exp=i, mod=17) for i in range(8)]), 8, 5, ValueError, _ERR_MUST_BE_CORRECT_ROOT),
]
TEST_CASES_Q17D8 = [
    (17, 8, 17//2, (17).bit_length())  # 0 modulus, 1 degree, 2 halfmod, 3 logmod
]
TEST_CASES_Q17D8 = [
    t + tuple([_find_primitive_root(modulus=t[0], root_order=2 * t[1])])  # 4 root
    for t in TEST_CASES_Q17D8
]
TEST_CASES_Q17D8 = [
    t + tuple([_bit_reverse_copy([pow(base=t[-1], exp=i, mod=t[0]) for i in range(t[1])])])  # 5 brv_root_powers
    for t in TEST_CASES_Q17D8
]
TEST_CASES_Q17D8 = [
    t + tuple([pow(base=t[4], exp=t[0]-2, mod=t[0])])  # 6 inv_root
    for t in TEST_CASES_Q17D8
]
TEST_CASES_Q17D8 = [
    t + tuple([_bit_reverse_copy([pow(base=t[-1], exp=i, mod=t[0]) for i in range(t[1])])])  # 7 brv_inv_root_powers
    for t in TEST_CASES_Q17D8
]
TEST_MONOMIALS = [
    [1, 0, 0, 0, 0, 0, 0, 0,], [0, 1, 0, 0, 0, 0, 0, 0,], [0, 0, 1, 0, 0, 0, 0, 0,], [0, 0, 0, 1, 0, 0, 0, 0,],
    [0, 0, 0, 0, 1, 0, 0, 0,], [0, 0, 0, 0, 0, 1, 0, 0,], [0, 0, 0, 0, 0, 0, 1, 0,], [0, 0, 0, 0, 0, 0, 0, 1,]]
TEST_FORWARD_MONOMIALS = [
    [1, 1, 1, 1, 1, 1, 1, 1], [3, -3, 5, -5, -7, 7, -6, 6], [-8, -8, 8, 8, -2, -2, 2, 2], [-7, 7, 6, -6, -3, 3, 5, -5],
    [-4, -4, -4, -4, 4, 4, 4, 4], [5, -5, -3, 3, 6, -6, -7, 7], [-2, -2, 2, 2, -8, -8, 8, 8], [-6, 6, -7, 7, 5, -5, 3, -3]]  # NTT(monomials)
TEST_BACKWARD_MONOMIALS = [
    [-2, 5, -4, -7, -8, 3, 1, 6], [-2, -5, -4, 7, -8, -3, 1, -6], [-2, 3, 4, -6, -8, -5, -1, -7], [-2, -3, 4, 6, -8, 5, -1, 7],
    [-2, -7, 1, -5, 8, -6, -4, 3], [-2, 7, 1, 5, 8, 6, -4, -3], [-2, 6, -1, 3, 8, -7, 4, 5], [-2, -6, -1, -3, 8, 7, 4, -5]]  # INTT(monomials)
TEST_MONOMIAL_MULTIPLICATION = [
    t + tuple([
        TEST_MONOMIALS[i],  # a unit vector
        TEST_FORWARD_MONOMIALS[i],  # the hand-computed forward transform of a
        TEST_BACKWARD_MONOMIALS[i], # the hand-computed inverse transform of a
    ])
    for t in TEST_CASES_Q17D8
    for i in range(len(TEST_MONOMIALS))
]
TEST_MONOMIAL_MULTIPLICATION = [
    t + tuple([
        _cooley_tukey_ntt(val=deepcopy(t[-3]), modulus=t[0], halfmod=t[2], logmod=t[3], root_order=2 * t[1],
                          bit_rev_root_powers=t[5]),  # forward transform of a
        _gentleman_sande_intt(val=deepcopy(t[-3]), modulus=t[0], halfmod=t[2], logmod=t[3], root_order=2 * t[1],
                              bit_rev_inv_root_powers=t[7]),  # inverse transform of a
        _gentleman_sande_intt(val=deepcopy(t[-2]), modulus=t[0], halfmod=t[2], logmod=t[3], root_order=2 * t[1],
                              bit_rev_inv_root_powers=t[7]),  # inverse transform of hand-computed forward transform
        _cooley_tukey_ntt(val=deepcopy(t[-1]), modulus=t[0], halfmod=t[2], logmod=t[3], root_order=2 * t[1],
                          bit_rev_root_powers=t[5]),   # forward transform of hand-computed inverse transform
    ])
    for t in TEST_MONOMIAL_MULTIPLICATION
]
TEST_RAN_MUL_TWO_RANDOM_POLYS = deepcopy(TEST_CASES_Q17D8)
TEST_RAN_MUL_TWO_RANDOM_POLYS = [
    t + tuple([
        [randbelow(t[0]) for _ in range(t[1])],
        [randbelow(t[0]) for _ in range(t[1])]
    ])
    for t in TEST_RAN_MUL_TWO_RANDOM_POLYS
    for _ in range(SAMPLE_SIZE)
]
TEST_RAN_MUL_CASES = []
for t in TEST_RAN_MUL_TWO_RANDOM_POLYS:
    f = t[-2]
    g = t[-1]
    fg = [0 for _ in range(2*t[1])]
    for i, a in enumerate(f):
        for j, b in enumerate(g):
            fg[i+j] += a*b
    fg = [_cent(val=x - y, modulus=t[0], halfmod=t[2], logmod=t[3]) for x, y in zip(fg[:t[1]], fg[t[1]:])]
    TEST_RAN_MUL_CASES += [t + tuple([
        fg,
        _cooley_tukey_ntt(val=deepcopy(f), modulus=t[0], halfmod=t[2], logmod=t[3], root_order=2 * t[1],
                          bit_rev_root_powers=t[5]),
        _cooley_tukey_ntt(val=deepcopy(g), modulus=t[0], halfmod=t[2], logmod=t[3], root_order=2 * t[1],
                          bit_rev_root_powers=t[5]),
        _cooley_tukey_ntt(val=deepcopy(fg), modulus=t[0], halfmod=t[2], logmod=t[3], root_order=2 * t[1],
                          bit_rev_root_powers=t[5]),
    ])]
TEST_RAN_MUL_CASES = [t + tuple([
    [_cent(val=x * y, modulus=t[0], halfmod=t[2], logmod=t[3]) for x, y in zip(t[-3], t[-2])]
]) for t in TEST_RAN_MUL_CASES]
TEST_RAN_MUL_CASES = [t + tuple([
    _gentleman_sande_intt(val=deepcopy(t[-1]), modulus=t[0], halfmod=t[2], logmod=t[3], root_order=2 * t[1],
                          bit_rev_inv_root_powers=t[7]),
]) for t in TEST_RAN_MUL_CASES]
TEST_RAN_MUL_CASES = [t + tuple([
    [_cent(val=x - y, modulus=t[0], halfmod=t[2], logmod=t[3]) for x, y in zip(t[-1], t[-6])]
]) for t in TEST_RAN_MUL_CASES]
TEST_RAN_MUL_CASES = [t + tuple([
    _ntt_poly_mult(f=t[8], g=t[9], modulus=t[0], halfmod=t[2], logmod=t[3], degree=t[1], root_order=2 * t[1], root=t[4], inv_root=t[6], brv_root_powers=t[5])
]) for t in TEST_RAN_MUL_CASES]


@pytest.mark.parametrize('number,expected_result', IS_ODD_PRIME_TEST_DATA)
def test_is_odd_prime(number, expected_result):
    assert _is_odd_prime(val=number) == expected_result


@pytest.mark.parametrize('modulus, root_order, expected_result', HAS_PRIMITIVE_ROOT_OF_UNITY_TEST_DATA)
def test_has_primitive_root_of_unity(modulus, root_order, expected_result):
    assert _has_primitive_root_of_unity(modulus=modulus, root_order=root_order) == expected_result


@pytest.mark.parametrize('number, expected_result', IS_POW_TWO_GEQ_TWO_TEST_DATA)
def test_is_pow_two_geq_two(number, expected_result):
    assert _is_pow_two_geq_two(val=number) == expected_result


@pytest.mark.parametrize('some_list, expected_result', BIT_REVERSE_COPY_TEST_DATA)
def test_bit_reverse_copy(some_list, expected_result):
    assert _bit_reverse_copy(val=some_list) == expected_result
    assert _bit_reverse_copy(val=expected_result) == some_list


@pytest.mark.parametrize("modulus,halfmod,logmod", CHECK_MODULUS_HALFMOD_LOGMOD_TEST_DATA)
def test_check_modulus_halfmod_logmod_valid(modulus, halfmod, logmod):
    _check_modulus_halfmod_logmod(modulus, halfmod, logmod)


@pytest.mark.parametrize("modulus,halfmod,logmod,exception_type,error_message", CHECK_MODULUS_HALFMOD_LOGMOD_ERRORS)
def test_check_modulus_halfmod_logmod_invalid(modulus, halfmod, logmod, exception_type, error_message):
    with pytest.raises(exception_type) as exc_info:
        _check_modulus_halfmod_logmod(modulus, halfmod, logmod)
    assert str(exc_info.value) == error_message


@pytest.mark.parametrize('val, modulus, halfmod, logmod, expected_value', CENT_TEST_DATA)
def test_cent(val, modulus, halfmod, logmod, expected_value):
    assert halfmod == modulus//2
    assert 2**(logmod-1) < modulus <= 2**logmod
    assert _cent(val=val, modulus=modulus, halfmod=halfmod, logmod=logmod) == expected_value


@pytest.mark.parametrize('val, modulus, root_order, expected_value', IS_ROOT_OF_UNITY_TEST_DATA)
def test_is_root_of_unity(val, modulus, root_order, expected_value):
    assert _is_root_of_unity(val=val, modulus=modulus, root_order=root_order) == expected_value


@pytest.mark.parametrize('val, modulus, root_order, expected_value', IS_PRIMITIVE_ROOT_OF_UNITY_TEST_DATA)
def test_is_primitive_root_of_unity(val, modulus, root_order, expected_value):
    assert _is_primitive_root(val=val, modulus=modulus, root_order=root_order) == expected_value


@pytest.mark.parametrize('modulus, root_order, expected_value', FIND_PRIMITIVE_ROOT_TEST_DATA)
def test_find_primitive_root(modulus, root_order, expected_value):
    assert _find_primitive_root(modulus=modulus, root_order=root_order) == expected_value


@pytest.mark.parametrize('val,modulus,root_order,inv_flag,root,brv_powers,halfmod,logmod', CHECK_NTT_AND_INTT_VALID_DATA)
def test_check_ntt_and_intt_valid_data(val, modulus, root_order, inv_flag, root, brv_powers, halfmod, logmod):
    _check_ntt_and_intt(val=val, modulus=modulus, halfmod=halfmod, logmod=logmod, root_order=root_order,
                        brv_powers=brv_powers, inv_flag=inv_flag)


@pytest.mark.parametrize('val,modulus,root_order,inv_flag,brv_powers,halfmod,logmod,exception_type,error_message', CHECK_NTT_AND_INTT_ERRORS)
def test_check_ntt_and_intt_errors(val, modulus, root_order, inv_flag, brv_powers, halfmod, logmod, exception_type, error_message):
    with pytest.raises(exception_type) as exc_info:
        _check_ntt_and_intt(val=val, modulus=modulus, halfmod=halfmod, logmod=logmod, root_order=root_order,
                            brv_powers=brv_powers, inv_flag=inv_flag)
    assert str(exc_info.value) == error_message


@pytest.mark.parametrize('modulus,degree,halfmod,logmod,root,brv_root_powers,inv_root,brv_inv_root_powers,' +
                         'input_vector,hand_ntt_input_vector,hand_intt_input_vector,ntt_input_vector,' +
                         'intt_input_vector,intt_hand_ntt_input_vector,ntt_hand_intt_input_vector', TEST_MONOMIAL_MULTIPLICATION)
def test_cooley_tukey_ntt_and_gentleman_sande_intt(modulus, degree, halfmod, logmod, root, brv_root_powers, inv_root, brv_inv_root_powers,
                                                   input_vector, hand_ntt_input_vector, hand_intt_input_vector, ntt_input_vector,
                                                   intt_input_vector, intt_hand_ntt_input_vector, ntt_hand_intt_input_vector):
    assert all((tomato - carrot) % modulus == 0 for tomato, carrot in zip(hand_ntt_input_vector, ntt_input_vector))
    assert all((tomato - carrot) % modulus == 0 for tomato, carrot in zip(hand_intt_input_vector, intt_input_vector))
    assert all((tomato - carrot) % modulus == 0 for tomato, carrot in zip(input_vector, intt_hand_ntt_input_vector))
    assert all((tomato - carrot) % modulus == 0 for tomato, carrot in zip(input_vector, ntt_hand_intt_input_vector))


@pytest.mark.parametrize('modulus,degree,halfmod,logmod,root,brv_root_powers,inv_root,brv_inv_root_powers,' +
                         'f,g,hand_fg,ntt_f,ntt_g,ntt_hand_fg,hadamard_fhat_ghat,intt_hadamard,diff,fg_polymul', TEST_RAN_MUL_CASES)
def test_ran_mul(modulus,degree,halfmod,logmod,root,brv_root_powers,inv_root,brv_inv_root_powers,f,g,hand_fg,ntt_f,ntt_g,ntt_hand_fg,hadamard_fhat_ghat,intt_hadamard,diff,fg_polymul):
    assert all(x % modulus == 0 for x in diff)
    assert all((tomato - carrot) % modulus == 0 for tomato, carrot in zip(hadamard_fhat_ghat,ntt_hand_fg))  # NTT(f) star NTT(g) = NTT(f * g)
    assert all((tomato - carrot) % modulus == 0 for tomato, carrot in zip(intt_hadamard, hand_fg))  # INTT(NTT(f) star NTT(g)) = f * g
    assert all((tomato - carrot) % modulus == 0 for tomato, carrot in zip(fg_polymul,hand_fg))


@pytest.mark.parametrize('modulus,degree,halfmod,logmod,root,brv_root_powers,inv_root,brv_inv_root_powers,' +
                         'f,g,hand_fg,ntt_f,ntt_g,ntt_hand_fg,hadamard_fhat_ghat,intt_hadamard,diff,fg_polymul', TEST_RAN_MUL_CASES)
def test_ntt_poly_mult(modulus,degree,halfmod,logmod,root,brv_root_powers,inv_root,brv_inv_root_powers,f,g,hand_fg,ntt_f,ntt_g,ntt_hand_fg,hadamard_fhat_ghat,intt_hadamard,diff,fg_polymul):
    observed_val = _ntt_poly_mult(f=f, g=g, modulus=modulus, halfmod=halfmod, logmod=logmod,
                                  degree=degree, root_order=2*degree, root=root, inv_root=inv_root,
                                  brv_root_powers=brv_root_powers)
    assert all((a-b) % modulus == 0 for a, b in zip(observed_val, hand_fg))
    assert all((a-b) % modulus == 0 for a, b in zip(observed_val, fg_polymul))


def test_derived_params():
    # simple, lazy
    pass
