import pytest
from copy import deepcopy
from random import randrange
from typing import List, Tuple
from algebra.ntt import find_prou, brv_copy, ntt_poly_mult
from algebra.errors import _INVALID_REP_TYPE_ERR
from algebra.polynomials import (
    _PolynomialRepresentation as PolyR,
    _PolynomialCoefficientRepresentation as PolyCR,
    _PolynomialNTTRepresentation as PolyNTTR,
    _POLYNOMIAL_COEFFICIENT_REPRESENTATION_STR_PREFIX,
    _POLYNOMIAL_NTT_REPRESENTATION_STR_PREFIX,
    Polynomial
)
from test_ntt import SAMPLE_SIZE, PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE

ARITHMETIC_TEST_CASES = [
    t + tuple([
        find_prou(mod=t[1], deg=t[0]),
    ])
    for t in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE]

ARITHMETIC_TEST_CASES = [
    t + tuple([
        [pow(base=t[-1], exp=i, mod=t[1]) for i in range(t[0])],
    ])
    for t in ARITHMETIC_TEST_CASES
]
ARITHMETIC_TEST_CASES = [
    t + tuple([
        brv_copy(t[-1]),
        pow(t[-2], t[1] - 2, t[1]),
    ])
    for t in ARITHMETIC_TEST_CASES
]
ARITHMETIC_TEST_CASES = [
    t + tuple([
        [pow(base=t[-1], exp=i, mod=t[1]) for i in range(t[0])]
    ])
    for t in ARITHMETIC_TEST_CASES
]
ARITHMETIC_TEST_CASES = [
    t + tuple([
        brv_copy(t[-1]),
        t[1] // 2,
        t[1].bit_length(),
    ])
    for t in ARITHMETIC_TEST_CASES
]


@pytest.mark.parametrize("deg,mod,root,root_powers,brv_powers,inv_root,inv_root_powers,brv_inv_root_powers,halfmod,logmod", ARITHMETIC_TEST_CASES)
def test_arithmetic(deg,mod,root,root_powers,brv_powers,inv_root,inv_root_powers,brv_inv_root_powers,halfmod,logmod):
    for _ in range(SAMPLE_SIZE):
        a_coefs: List[int] = [randrange(mod) for _ in range(deg)]
        b_coefs: List[int] = [randrange(mod) for _ in range(deg)]
        a: PolyCR = PolyCR(
            mod=mod,
            vals=tuple(a_coefs),
        )
        b: PolyCR = PolyCR(
            mod=mod,
            vals=tuple(b_coefs),
        )
        expected_c_coefs: List[int] = [0 for _ in range(2 * deg)]
        # foil
        for i, next_a in enumerate(a_coefs):
            for j, next_b in enumerate(b_coefs):
                expected_c_coefs[i + j] = (
                    expected_c_coefs[i + j] + next_a * next_b
                ) % mod
        expected_c_coefs = [
            (x - y) % mod for x, y in zip(expected_c_coefs[:deg], expected_c_coefs[deg:])
        ]
        expected_c: PolyCR = PolyCR(
            mod=mod,
            vals=tuple(expected_c_coefs),
        )

        observed_c: PolyCR = a * b
        assert len(observed_c._vals) == len(expected_c._vals)
        assert all(
            (x - y) % mod == 0
            for x, y in zip(observed_c._vals, expected_c._vals)
        )
        assert expected_c == observed_c
        dc_a_coefs: List[int] = deepcopy(a_coefs)
        dc_b_coefs: List[int] = deepcopy(b_coefs)
        another_observed_c_coefs: List[int] = ntt_poly_mult(f=dc_a_coefs, g=dc_b_coefs, mod=mod)
        right: List[int] = [-y for y in another_observed_c_coefs[deg:]]
        another_observed_c_coefs = another_observed_c_coefs[:deg]
        while len(right) >= deg:
            next_right = right[:deg]
            right = [-y for y in right[deg:]]
            for i, y in enumerate(next_right):
                another_observed_c_coefs[i] = (another_observed_c_coefs[i] + y) % mod
        for i, y in enumerate(right):
            another_observed_c_coefs[i] = (another_observed_c_coefs[i] + y) % mod
        assert len(another_observed_c_coefs) == len(expected_c._vals)
        assert all(
            (x - y) % mod == 0
            for x, y in zip(observed_c._vals, another_observed_c_coefs))
        assert all(
            (x - y) % mod == 0
            for x, y in zip(expected_c._vals, observed_c._vals))
        another_observed_c: PolyCR = PolyCR(mod=mod, vals=tuple(another_observed_c_coefs))
        assert another_observed_c == observed_c
        # Since expected_c and observed_c are equivalent, we only need one in the rest of this test.

        a_hat: PolyNTTR = a.transform_to_ntt_rep()
        b_hat: PolyNTTR = b.transform_to_ntt_rep()
        c_hat: PolyNTTR = observed_c.transform_to_ntt_rep()
        a_hat_times_b_hat: PolyNTTR = a_hat * b_hat
        assert c_hat == a_hat_times_b_hat

        transformed_c_hat: PolyCR = c_hat.transform_to_coef_rep()
        transformed_a_hat_times_b_hat: PolyCR = a_hat_times_b_hat.transform_to_coef_rep()
        assert observed_c == transformed_c_hat
        assert observed_c == transformed_a_hat_times_b_hat

        inv_a_hat: PolyCR = a_hat.transform_to_coef_rep()
        inv_b_hat: PolyCR = b_hat.transform_to_coef_rep()
        inv_a_hat_times_b_hat: PolyCR = a_hat_times_b_hat.transform_to_coef_rep()
        inv_c_hat: PolyCR = c_hat.transform_to_coef_rep()
        assert inv_a_hat == a
        assert inv_b_hat == b
        assert inv_a_hat_times_b_hat == observed_c
        assert inv_a_hat_times_b_hat == inv_c_hat


def test_monomial_products():
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        root: int = find_prou(mod=q, deg=d)
        inv_root: int = pow(root, q - 2, q)
        for _ in range(SAMPLE_SIZE):
            a_coef: int = randrange(1, q)
            a_index: int = randrange(d)
            a_coefs: List[int] = [0 if i != a_index else a_coef for i in range(d)]
            a: PolyCR = PolyCR(
                mod=q,
                vals=tuple(a_coefs))

            b_coef: int = randrange(1, q)
            b_index: int = randrange(d)
            b_coefs: List[int] = [0 if i != b_index else b_coef for i in range(d)]
            b: PolyCR = PolyCR(
                mod=q,
                vals=tuple(b_coefs))

            expected_c_coefs: List[int] = [0 for _ in range(d)]
            expected_c_coefs[(a_index + b_index) % d] = ((a_coef * b_coef) % q) * (
                1 - 2 * int(a_index + b_index >= d)
            )
            expected_c: PolyCR = PolyCR(
                mod=q,
                vals=tuple(expected_c_coefs))
            observed_c: PolyCR = a * b
            try:
                assert expected_c == observed_c
            except AssertionError:
                expected_c: PolyCR = PolyCR(
                    mod=q,
                    vals=tuple(expected_c_coefs),
                )
                observed_c: PolyCR = a * b
                assert expected_c == observed_c

            a_hat: PolyNTTR = a.transform_to_ntt_rep()
            b_hat: PolyNTTR = b.transform_to_ntt_rep()
            a_hat_times_b_hat: PolyNTTR = a_hat * b_hat
            inv_a_hat_times_b_hat: PolyCR = a_hat_times_b_hat.transform_to_coef_rep()
            assert inv_a_hat_times_b_hat == observed_c


def test_poly_init():
    with pytest.raises(TypeError):
        PolyCR(mod=1, vals=1)
    with pytest.raises(TypeError):
        PolyCR(mod=5, vals=1)
    with pytest.raises(TypeError):
        PolyCR(mod=1.0, vals=["hello world"])
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        root: int = find_prou(mod=q, deg=d)
        inv_root: int = pow(root, q - 2, q)
        for _ in range(SAMPLE_SIZE):
            a_coefs: Tuple[int, ...] = tuple([randrange(q) for _ in range(d)])
            a: PolyCR = PolyCR(
                mod=q,
                vals=tuple(a_coefs))
            assert a.modulus == q
            assert a.deg == d
            assert a.root_order == 2 * d
            assert a.root == root
            assert a.inv_root == inv_root
            assert a._vals == a_coefs


def test_poly_str():
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        root: int = find_prou(mod=q, deg=d)
        inv_root: int = pow(root, q - 2, q)
        for _ in range(SAMPLE_SIZE):
            a_coefs: Tuple[int, ...] = tuple([randrange(q) for _ in range(d)])
            a: PolyCR = PolyCR(
                mod=q,
                vals=tuple(a_coefs))
            assert (
                str(a)
                == _POLYNOMIAL_COEFFICIENT_REPRESENTATION_STR_PREFIX + f"(modulus={q}, degree={d}, root={root}, inv_root={inv_root}, values={a_coefs})")


def test_poly_eq():
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        for _ in range(SAMPLE_SIZE):
            a_coefs: List[int] = [randrange(q) for _ in range(d)]
            a: PolyCR = PolyCR(
                mod=q,
                vals=tuple(a_coefs))
            b: PolyCR = PolyCR(
                mod=q,
                vals=tuple(a._vals))
            assert a == a
            assert a == b
            assert len(a._vals) == len(b._vals)


def test_poly_add():
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        root: int = find_prou(mod=q, deg=d)
        inv_root: int = pow(root, q - 2, q)
        for _ in range(SAMPLE_SIZE):
            a_coefs: List[int] = [randrange(1, q) for _ in range(d)]
            b_coefs: List[int] = [randrange(1, q) for _ in range(d)]
            a: PolyCR = PolyCR(
                mod=q,
                vals=tuple(a_coefs))
            b: PolyCR = PolyCR(
                mod=q,
                vals=tuple(b_coefs))
            c: PolyCR = a + b
            assert c.modulus == q
            assert c.deg == d
            assert c.root_order == 2 * d
            assert c.root == root
            assert c.inv_root == inv_root
            for i in range(d):
                assert (c._vals[i] - (a._vals[i] + b._vals[i])) % q == 0


def test_poly_sub():
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        root: int = find_prou(mod=q, deg=d)
        inv_root: int = pow(base=root, exp=q-2, mod=q)
        for _ in range(SAMPLE_SIZE):
            a_coefs: List[int] = [randrange(1, q) for _ in range(d)]
            b_coefs: List[int] = [randrange(1, q) for _ in range(d)]
            a: PolyCR = PolyCR(
                mod=q,
                vals=tuple(a_coefs))
            b: PolyCR = PolyCR(
                mod=q,
                vals=tuple(b_coefs))
            c: PolyCR = a - b
            assert c.modulus == q
            assert c.deg == d
            assert c.root_order == 2 * d
            assert c.root == root
            assert c.inv_root == inv_root
            for i in range(d):
                assert (c._vals[i] - (a._vals[i] - b._vals[i])) % q == 0


def test_poly_mul():
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        root: int = find_prou(mod=q, deg=d)
        inv_root: int = pow(root, q - 2, q)
        for _ in range(SAMPLE_SIZE):
            a_coefs: List[int] = [randrange(1, q) for _ in range(d)]
            b_coefs: List[int] = [randrange(1, q) for _ in range(d)]
            expected_c_coefs: List[int] = [0 for _ in range(2 * d)]
            for i, next_a_coef in enumerate(a_coefs):
                for j, next_b_coef in enumerate(b_coefs):
                    expected_c_coefs[i + j] = (
                        expected_c_coefs[i + j] + (next_a_coef * next_b_coef)
                    ) % q
            expected_c_coefs = [
                (x - y) % q for x, y in zip(expected_c_coefs[:d], expected_c_coefs[d:])
            ]

            a: PolyCR = PolyCR(
                mod=q,
                vals=tuple(a_coefs))
            b: PolyCR = PolyCR(
                mod=q,
                vals=tuple(b_coefs))
            expected_c: PolyCR = PolyCR(
                mod=q,
                vals=tuple(expected_c_coefs))
            observed_c: PolyCR = a * b
            assert expected_c == observed_c


def test_poly_norm():
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        for _ in range(SAMPLE_SIZE):
            a_coefs: List[int] = [randrange(1, q) for _ in range(d)]
            a: PolyCR = PolyCR(
                mod=q,
                vals=tuple(a_coefs))
            expected_infinity_norm: int = max(abs(x) for x in a._vals)
            observed_infinity_norm_and_weight: Tuple[int, int] = a.coef_norm_wght
            assert expected_infinity_norm, expected_infinity_norm == observed_infinity_norm_and_weight


def test_poly_ntt_init():
    with pytest.raises(TypeError):
        PolyNTTR(mod=2, vals=1)
    with pytest.raises(TypeError):
        PolyNTTR(mod=5, vals=1)
    with pytest.raises(TypeError):
        PolyNTTR(mod=5, vals=["hello world"])
    with pytest.raises(TypeError):
        PolyNTTR(mod=5, vals=["hello world"])
    assert PolyNTTR(mod=5, vals=tuple([1]))
    with pytest.raises(TypeError):
        PolyNTTR(mod=5, vals=[1, 2, 3])
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        for _ in range(SAMPLE_SIZE):
            a_coefs: List[int] = [randrange(1, q) for _ in range(d)]
            a: PolyCR = PolyCR(
                mod=q,
                vals=tuple(a_coefs))
            b: PolyCR = PolyCR(
                mod=q,
                vals=tuple(a._vals))
            assert a == b
            assert a == a


def test_poly_ntt_str():
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        root: int = find_prou(mod=q, deg=d)
        inv_root: int = pow(root, q - 2, q)
        for _ in range(SAMPLE_SIZE):
            a_vals: List[int] = [randrange(1, q) for _ in range(d)]
            a_hat: PolyNTTR = PolyNTTR(
                mod=q,
                vals=tuple(a_vals))
            assert (
                str(a_hat)
                == _POLYNOMIAL_NTT_REPRESENTATION_STR_PREFIX+f"(modulus={q}, degree={d}, root={root}, inv_root={inv_root}, values={tuple(a_vals)})")


def test_poly_ntt_eq():
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        root: int = find_prou(mod=q, deg=d)
        inv_root: int = pow(root, q - 2, q)
        for _ in range(SAMPLE_SIZE):
            a_vals: List[int] = [randrange(1, q) for _ in range(d)]
            a_hat: PolyNTTR = PolyNTTR(
                mod=q,
                vals=tuple(a_vals))
            b_hat: PolyNTTR = PolyNTTR(
                mod=q,
                vals=tuple(a_hat._vals))
            assert a_hat == b_hat
            assert a_hat == a_hat


def test_poly_ntt_add():
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        root: int = find_prou(mod=q, deg=d)
        inv_root: int = pow(root, q - 2, q)
        for _ in range(SAMPLE_SIZE):
            a_vals: List[int] = [randrange(1, q) for _ in range(d)]
            a_hat: PolyNTTR = PolyNTTR(
                mod=q,
                vals=tuple(a_vals))
            b_vals: List[int] = [randrange(1, q) for _ in range(d)]
            b_hat: PolyNTTR = PolyNTTR(
                mod=q,
                vals=tuple(b_vals))
            c_hat: PolyNTTR = a_hat + b_hat
            assert all([
                (z - (x + y)) % q == 0 for x, y, z in zip(a_vals, b_vals, c_hat._vals)
            ])
            with pytest.raises(TypeError):
                c_hat: PolyNTTR = a_hat + b_vals
            with pytest.raises(TypeError):
                c_hat: PolyNTTR = a_vals + b_hat
            with pytest.raises(TypeError):
                c_hat: PolyNTTR = b_hat + a_vals
            with pytest.raises(TypeError):
                c_hat: PolyNTTR = a_hat + b_vals[0]
            with pytest.raises(TypeError):
                c_hat: PolyNTTR = a_vals[0] + b_hat
            with pytest.raises(TypeError):
                c_hat: PolyNTTR = b_hat + a_vals[0]
            with pytest.raises(TypeError):
                c_hat: PolyNTTR = a_hat + 0
            c_hat: PolyNTTR = 0 + a_hat
            with pytest.raises(TypeError):
                c_hat: PolyNTTR = a_hat + 1


def test_poly_ntt_sub():
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        root: int = find_prou(mod=q, deg=d)
        inv_root: int = pow(root, q - 2, q)
        for _ in range(SAMPLE_SIZE):
            a_vals: List[int] = [randrange(1, q) for _ in range(d)]
            a_hat: PolyNTTR = PolyNTTR(
                mod=q,
                vals=tuple(a_vals))
            b_vals: List[int] = [randrange(1, q) for _ in range(d)]
            b_hat: PolyNTTR = PolyNTTR(
                mod=q,
                vals=tuple(b_vals))
            c_hat: PolyNTTR = a_hat - b_hat
            for i in range(d):
                assert (c_hat._vals[i] - (a_vals[i] - b_vals[i])) % q == 0
            with pytest.raises(TypeError):
                c_hat: PolyNTTR = a_hat - b_vals
            with pytest.raises(TypeError):
                c_hat: PolyNTTR = a_vals - b_hat
            with pytest.raises(TypeError):
                c_hat: PolyNTTR = a_hat - b_vals
            with pytest.raises(TypeError):
                c_hat: PolyNTTR = a_hat - b_vals[0]
            with pytest.raises(TypeError):
                c_hat: PolyNTTR = a_vals[0] - b_hat
            with pytest.raises(TypeError):
                c_hat: PolyNTTR = a_hat - 0
                assert c_hat == a_hat
            with pytest.raises(TypeError):
                c_hat: PolyNTTR = 0 - a_hat
                assert c_hat == -a_hat
            with pytest.raises(TypeError):
                c_hat: PolyNTTR = a_hat - 1


def test_poly_ntt_neg():
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        root: int = find_prou(mod=q, deg=d)
        inv_root: int = pow(root, q - 2, q)
        for _ in range(SAMPLE_SIZE):
            a_vals: List[int] = [randrange(1, q) for _ in range(d)]
            a_hat: PolyNTTR = PolyNTTR(
                mod=q,
                vals=tuple(a_vals))
            b_hat: PolyNTTR = -a_hat
            for i in range(d):
                assert (b_hat._vals[i] + a_vals[i]) % q == 0


def test_poly_ntt_radd():
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        root: int = find_prou(mod=q, deg=d)
        inv_root: int = pow(root, q - 2, q)
        for _ in range(SAMPLE_SIZE):
            a_vals: List[int] = [randrange(1, q) for _ in range(d)]
            a_hat: PolyNTTR = PolyNTTR(
                mod=q,
                vals=tuple(a_vals))
            b_vals: List[int] = [randrange(1, q) for _ in range(d)]
            b_hat: PolyNTTR = PolyNTTR(
                mod=q,
                vals=tuple(b_vals))
            c_hat: PolyNTTR = b_hat + a_hat
            for i in range(d):
                assert (c_hat._vals[i] - (a_vals[i] + b_vals[i])) % q == 0
            with pytest.raises(TypeError):
                c_hat: PolyNTTR = b_hat + a_vals
            with pytest.raises(TypeError):
                c_hat: PolyNTTR = b_vals + a_hat
            with pytest.raises(TypeError):
                c_hat: PolyNTTR = b_hat + a_vals[0]
            with pytest.raises(TypeError):
                c_hat: PolyNTTR = b_vals[0] + a_hat
            c_hat: PolyNTTR = 0 + b_hat
            assert c_hat == b_hat
            with pytest.raises(TypeError):
                c_hat: PolyNTTR = b_hat + q


def test_poly_ntt_mul():
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        root: int = find_prou(mod=q, deg=d)
        inv_root: int = pow(root, q - 2, q)
        for _ in range(SAMPLE_SIZE):
            a_vals: List[int] = [randrange(1, q) for _ in range(d)]
            a_hat: PolyNTTR = PolyNTTR(
                mod=q,
                vals=tuple(a_vals))
            b_vals: List[int] = [randrange(1, q) for _ in range(d)]
            b_hat: PolyNTTR = PolyNTTR(
                mod=q,
                vals=tuple(b_vals))
            c_hat: PolyNTTR = a_hat * b_hat
            for i in range(d):
                assert (c_hat._vals[i] - (a_vals[i] * b_vals[i])) % q == 0
            with pytest.raises(TypeError):
                a_vals * b_hat
            with pytest.raises(TypeError):
                a_hat * b_vals
            with pytest.raises(TypeError):
                c_hat: PolyNTTR = a_hat * 0
            with pytest.raises(TypeError):
                c_hat: int = 0 * b_hat
            with pytest.raises(TypeError):
                c_hat: PolyNTTR = a_hat * 1


def test_poly_ntt_rmul():
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        root: int = find_prou(mod=q, deg=d)
        inv_root: int = pow(root, q - 2, q)
        for _ in range(SAMPLE_SIZE):
            a_vals: List[int] = [randrange(1, q) for _ in range(d)]
            a_hat: PolyNTTR = PolyNTTR(
                mod=q,
                vals=tuple(a_vals))
            b_vals: List[int] = [randrange(1, q) for _ in range(d)]
            b_hat: PolyNTTR = PolyNTTR(
                mod=q,
                vals=tuple(b_vals))
            c_hat: PolyNTTR = b_hat * a_hat
            for i in range(d):
                assert (c_hat._vals[i] - (a_vals[i] * b_vals[i])) % q == 0
            with pytest.raises(TypeError):
                b_hat * a_vals
            with pytest.raises(TypeError):
                b_vals * a_hat
            with pytest.raises(TypeError):
                c_hat: PolyNTTR = b_hat * 0
            # assert c_hat == 0
            with pytest.raises(TypeError):
                c_hat: int = 0 * a_hat


def test_transform_2d():
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        root: int = find_prou(mod=q, deg=d)
        inv_root: int = pow(root, q - 2, q)
        for _ in range(SAMPLE_SIZE):
            a_coefs: List[int] = [randrange(1, q) for _ in range(d)]
            a: PolyCR = PolyCR(
                mod=q,
                vals=tuple(a_coefs))
            b_coefs: List[int] = [randrange(1, q) for _ in range(d)]
            b: PolyCR = PolyCR(
                mod=q,
                vals=tuple(b_coefs))
            c_coefs: List[int] = [0 for _ in range(2 * d)]
            for i in range(len(a_coefs)):
                for j in range(len(b_coefs)):
                    c_coefs[i + j] = (c_coefs[i + j] + a_coefs[i] * b_coefs[j]) % q
            c_coefs = [(x - y) % q for x, y in zip(c_coefs[:d], c_coefs[d:])]
            c: PolyCR = PolyCR(
                mod=q,
                vals=tuple(c_coefs))

            a_ntt: PolyNTTR = a.transform_to_ntt_rep()
            b_ntt: PolyNTTR = b.transform_to_ntt_rep()
            a_ntt_times_b_ntt: PolyNTTR = a_ntt * b_ntt
            inv_a_ntt_times_b_ntt: PolyCR = a_ntt_times_b_ntt.transform_to_coef_rep()
            assert inv_a_ntt_times_b_ntt == c


def test_comprehensive():
    for _ in range(SAMPLE_SIZE):
        for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
            root: int = find_prou(mod=q, deg=d)
            inv_root: int = pow(root, q - 2, q)
            # Random polynomial a
            a_coefs: List[int] = [randrange(1, q) for _ in range(d)]
            a: PolyCR = PolyCR(
                mod=q,
                vals=tuple(a_coefs))
            # Random polynomial b
            b_coefs: List[int] = [randrange(1, q) for _ in range(d)]
            b: PolyCR = PolyCR(
                mod=q,
                vals=tuple(b_coefs))

            # Transformed a and b
            a_hat: PolyNTTR = a.transform_to_ntt_rep()
            b_hat: PolyNTTR = b.transform_to_ntt_rep()

            # Product of transforms
            a_hat_times_b_hat: PolyNTTR = a_hat * b_hat

            # Invert
            inv_a_hat_times_b_hat: PolyCR = a_hat_times_b_hat.transform_to_coef_rep()

            # Check
            assert inv_a_hat_times_b_hat == a * b

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
    c,n,w = p.coef_norm_wght
    assert isinstance(p, Polynomial)
    assert p.modulus == MODULUS
    assert c == COEFFICIENTS  # Assuming transformation to NTT keeps values same in this example.


def test_initialization_with_ntt():
    p = Polynomial(MODULUS, COEFFICIENTS, NTT_REP)
    assert isinstance(p, Polynomial)
    assert p.modulus == MODULUS
    assert p.vals == COEFFICIENTS


def test_initialization_with_invalid_representation():
    with pytest.raises(ValueError) as e:
        Polynomial(MODULUS, COEFFICIENTS, INVALID_REP)
    assert str(e.value) == _INVALID_REP_TYPE_ERR


def test_polynomial_properties(poly_ntt):
    assert poly_ntt.modulus == MODULUS
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
    expected_representation = f"_PolynomialRepresentation(ntt={poly_ntt.ntt_rep})"
    assert str(poly_ntt) == expected_representation


@pytest.fixture
def example_polynomial_1() -> PolyR:
    """Example polynomial representation for testing."""
    mod = 17
    vals = (6, 3, 9, 10)
    return PolyR(mod=mod, vals=vals)


def test_to_bytes(example_polynomial_1: PolyR):
    """Test serialization of polynomial to bytes."""
    serialized_bytes = example_polynomial_1.to_bytes()
    expected_byte_length = 4 + 1 + 4 * example_polynomial_1.deg
    assert isinstance(serialized_bytes, bytes)
    assert len(serialized_bytes) == expected_byte_length
    mod, d = serialized_bytes[:5]
    assert mod == example_polynomial_1.modulus
    assert 1 << d == example_polynomial_1.deg


def test_from_bytes(example_polynomial_1: PolyR):
    """Test deserialization of polynomial from bytes."""
    serialized_bytes = example_polynomial_1.to_bytes()
    deserialized_poly = PolyR.from_bytes(data=serialized_bytes, representation_format=PolyR)
    assert deserialized_poly == example_polynomial_1
    assert deserialized_poly.modulus == example_polynomial_1.modulus
    assert deserialized_poly.vals == example_polynomial_1.vals


def test_to_from_bytes_round_trip(example_polynomial_1: PolyR):
    """Test serialization and deserialization is a lossless round-trip."""
    serialized_bytes = example_polynomial_1.to_bytes()
    deserialized_poly = PolyR.from_bytes(data=serialized_bytes, representation_format=PolyR)
    assert deserialized_poly == example_polynomial_1
    assert example_polynomial_1.to_bytes() == serialized_bytes
