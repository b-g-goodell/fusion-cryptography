import pytest
from copy import deepcopy
from random import randrange
from typing import List, Tuple
from api.ntt import find_prou, bit_reverse_copy, ntt_poly_mult
from algebra.polynomials import (
    PolynomialCoefficientRepresentation as Poly,
    PolynomialNTTRepresentation as PolyNTT,
    _POLYNOMIAL_COEFFICIENT_REPRESENTATION_STR_PREFIX,
    _POLYNOMIAL_NTT_REPRESENTATION_STR_PREFIX
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
        bit_reverse_copy(t[-1]),
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
        bit_reverse_copy(t[-1]),
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
        a: Poly = Poly(
            mod=mod,
            vals=tuple(a_coefs),
        )
        b: Poly = Poly(
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
        expected_c: Poly = Poly(
            mod=mod,
            vals=tuple(expected_c_coefs),
        )

        observed_c: Poly = a * b
        assert len(observed_c.vals) == len(expected_c.vals)
        assert all(
            (x - y) % mod == 0
            for x, y in zip(observed_c.vals, expected_c.vals)
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
        assert len(another_observed_c_coefs) == len(expected_c.vals)
        assert all(
            (x - y) % mod == 0
            for x, y in zip(observed_c.vals, another_observed_c_coefs))
        assert all(
            (x - y) % mod == 0
            for x, y in zip(expected_c.vals, observed_c.vals))
        another_observed_c: Poly = Poly(mod=mod, vals=tuple(another_observed_c_coefs))
        assert another_observed_c == observed_c
        # Since expected_c and observed_c are equivalent, we only need one in the rest of this test.

        a_hat: PolyNTT = a.transform()
        b_hat: PolyNTT = b.transform()
        c_hat: PolyNTT = observed_c.transform()
        a_hat_times_b_hat: PolyNTT = a_hat * b_hat
        assert c_hat == a_hat_times_b_hat

        transformed_c_hat: Poly = c_hat.transform()
        transformed_a_hat_times_b_hat: Poly = a_hat_times_b_hat.transform()
        assert observed_c == transformed_c_hat
        assert observed_c == transformed_a_hat_times_b_hat

        inv_a_hat: Poly = a_hat.transform()
        inv_b_hat: Poly = b_hat.transform()
        inv_a_hat_times_b_hat: Poly = a_hat_times_b_hat.transform()
        inv_c_hat: Poly = c_hat.transform()
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
            a: Poly = Poly(
                mod=q,
                vals=tuple(a_coefs))

            b_coef: int = randrange(1, q)
            b_index: int = randrange(d)
            b_coefs: List[int] = [0 if i != b_index else b_coef for i in range(d)]
            b: Poly = Poly(
                mod=q,
                vals=tuple(b_coefs))

            expected_c_coefs: List[int] = [0 for _ in range(d)]
            expected_c_coefs[(a_index + b_index) % d] = ((a_coef * b_coef) % q) * (
                1 - 2 * int(a_index + b_index >= d)
            )
            expected_c: Poly = Poly(
                mod=q,
                vals=tuple(expected_c_coefs))
            observed_c: Poly = a * b
            try:
                assert expected_c == observed_c
            except AssertionError:
                expected_c: Poly = Poly(
                    mod=q,
                    vals=tuple(expected_c_coefs),
                )
                observed_c: Poly = a * b
                assert expected_c == observed_c

            a_hat: PolyNTT = a.transform()
            b_hat: PolyNTT = b.transform()
            a_hat_times_b_hat: PolyNTT = a_hat * b_hat
            inv_a_hat_times_b_hat: Poly = a_hat_times_b_hat.transform()
            assert inv_a_hat_times_b_hat == observed_c


def test_poly_init():
    with pytest.raises(TypeError):
        Poly(mod=1, vals=1)
    with pytest.raises(TypeError):
        Poly(mod=5, vals=1)
    with pytest.raises(TypeError):
        Poly(mod=1.0, vals=["hello world"])
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        root: int = find_prou(mod=q, deg=d)
        inv_root: int = pow(root, q - 2, q)
        for _ in range(SAMPLE_SIZE):
            a_coefs: Tuple[int, ...] = tuple([randrange(q) for _ in range(d)])
            a: Poly = Poly(
                mod=q,
                vals=tuple(a_coefs))
            assert a.mod == q
            assert a.deg == d
            assert a.root_order == 2 * d
            assert a.root == root
            assert a.inv_root == inv_root
            assert a.vals == a_coefs


def test_poly_str():
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        root: int = find_prou(mod=q, deg=d)
        inv_root: int = pow(root, q - 2, q)
        for _ in range(SAMPLE_SIZE):
            a_coefs: Tuple[int, ...] = tuple([randrange(q) for _ in range(d)])
            a: Poly = Poly(
                mod=q,
                vals=tuple(a_coefs))
            assert (
                str(a)
                == _POLYNOMIAL_COEFFICIENT_REPRESENTATION_STR_PREFIX + f"(modulus={q}, degree={d}, root={root}, inv_root={inv_root}, values={a_coefs})")


def test_poly_repr():
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        root: int = find_prou(mod=q, deg=d)
        inv_root: int = pow(root, q - 2, q)
        for _ in range(SAMPLE_SIZE):
            a_coefs: List[int] = [randrange(q) for _ in range(d)]
            a: Poly = Poly(
                mod=q,
                vals=tuple(a_coefs))
            r = repr(a)
            s = _POLYNOMIAL_COEFFICIENT_REPRESENTATION_STR_PREFIX + f"(modulus={q}, degree={d}, root={root}, inv_root={inv_root}, values={tuple(a_coefs)})"
            assert r == s


def test_poly_eq():
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        for _ in range(SAMPLE_SIZE):
            a_coefs: List[int] = [randrange(q) for _ in range(d)]
            a: Poly = Poly(
                mod=q,
                vals=tuple(a_coefs))
            b: Poly = Poly(
                mod=q,
                vals=tuple(a.vals))
            assert a == a
            assert a == b
            assert len(a.vals) == len(b.vals)


def test_poly_add():
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        root: int = find_prou(mod=q, deg=d)
        inv_root: int = pow(root, q - 2, q)
        for _ in range(SAMPLE_SIZE):
            a_coefs: List[int] = [randrange(1, q) for _ in range(d)]
            b_coefs: List[int] = [randrange(1, q) for _ in range(d)]
            a: Poly = Poly(
                mod=q,
                vals=tuple(a_coefs))
            b: Poly = Poly(
                mod=q,
                vals=tuple(b_coefs))
            c: Poly = a + b
            assert c.mod == q
            assert c.deg == d
            assert c.root_order == 2 * d
            assert c.root == root
            assert c.inv_root == inv_root
            for i in range(d):
                assert (c.vals[i] - (a.vals[i] + b.vals[i])) % q == 0


def test_poly_sub():
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        root: int = find_prou(mod=q, deg=d)
        inv_root: int = pow(base=root, exp=q-2, mod=q)
        for _ in range(SAMPLE_SIZE):
            a_coefs: List[int] = [randrange(1, q) for _ in range(d)]
            b_coefs: List[int] = [randrange(1, q) for _ in range(d)]
            a: Poly = Poly(
                mod=q,
                vals=tuple(a_coefs))
            b: Poly = Poly(
                mod=q,
                vals=tuple(b_coefs))
            c: Poly = a - b
            assert c.mod == q
            assert c.deg == d
            assert c.root_order == 2 * d
            assert c.root == root
            assert c.inv_root == inv_root
            for i in range(d):
                assert (c.vals[i] - (a.vals[i] - b.vals[i])) % q == 0


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

            a: Poly = Poly(
                mod=q,
                vals=tuple(a_coefs))
            b: Poly = Poly(
                mod=q,
                vals=tuple(b_coefs))
            expected_c: Poly = Poly(
                mod=q,
                vals=tuple(expected_c_coefs))
            observed_c: Poly = a * b
            assert expected_c == observed_c


def test_poly_norm():
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        for _ in range(SAMPLE_SIZE):
            a_coefs: List[int] = [randrange(1, q) for _ in range(d)]
            a: Poly = Poly(
                mod=q,
                vals=tuple(a_coefs))
            expected_infinity_norm: int = max(abs(x) for x in a.vals)
            observed_infinity_norm_and_weight: Tuple[int, int] = a.coefs_norm_wght
            assert expected_infinity_norm, expected_infinity_norm == observed_infinity_norm_and_weight


def test_poly_ntt_init():
    with pytest.raises(TypeError):
        PolyNTT(mod=2, vals=1)
    with pytest.raises(TypeError):
        PolyNTT(mod=5, vals=1)
    with pytest.raises(TypeError):
        PolyNTT(mod=5, vals=["hello world"])
    with pytest.raises(TypeError):
        PolyNTT(mod=5, vals=["hello world"])
    assert PolyNTT(mod=5, vals=tuple([1]))
    with pytest.raises(TypeError):
        PolyNTT(mod=5, vals=[1, 2, 3])
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        for _ in range(SAMPLE_SIZE):
            a_coefs: List[int] = [randrange(1, q) for _ in range(d)]
            a: Poly = Poly(
                mod=q,
                vals=tuple(a_coefs))
            b: Poly = Poly(
                mod=q,
                vals=tuple(a.vals))
            assert a == b
            assert a == a


def test_poly_ntt_str():
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        root: int = find_prou(mod=q, deg=d)
        inv_root: int = pow(root, q - 2, q)
        for _ in range(SAMPLE_SIZE):
            a_vals: List[int] = [randrange(1, q) for _ in range(d)]
            a_hat: PolyNTT = PolyNTT(
                mod=q,
                vals=tuple(a_vals))
            assert (
                str(a_hat)
                == repr(a_hat)
                == _POLYNOMIAL_NTT_REPRESENTATION_STR_PREFIX+f"(modulus={q}, degree={d}, root={root}, inv_root={inv_root}, values={tuple(a_vals)})")


def test_poly_ntt_eq():
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        root: int = find_prou(mod=q, deg=d)
        inv_root: int = pow(root, q - 2, q)
        for _ in range(SAMPLE_SIZE):
            a_vals: List[int] = [randrange(1, q) for _ in range(d)]
            a_hat: PolyNTT = PolyNTT(
                mod=q,
                vals=tuple(a_vals))
            b_hat: PolyNTT = PolyNTT(
                mod=q,
                vals=tuple(a_hat.vals))
            assert a_hat == b_hat
            assert a_hat == a_hat


def test_poly_ntt_add():
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        root: int = find_prou(mod=q, deg=d)
        inv_root: int = pow(root, q - 2, q)
        for _ in range(SAMPLE_SIZE):
            a_vals: List[int] = [randrange(1, q) for _ in range(d)]
            a_hat: PolyNTT = PolyNTT(
                mod=q,
                vals=tuple(a_vals))
            b_vals: List[int] = [randrange(1, q) for _ in range(d)]
            b_hat: PolyNTT = PolyNTT(
                mod=q,
                vals=tuple(b_vals))
            c_hat: PolyNTT = a_hat + b_hat
            assert all([
                (z - (x + y)) % q == 0 for x, y, z in zip(a_vals, b_vals, c_hat.vals)
            ])
            with pytest.raises(TypeError):
                c_hat: PolyNTT = a_hat + b_vals
            with pytest.raises(TypeError):
                c_hat: PolyNTT = a_vals + b_hat
            with pytest.raises(TypeError):
                c_hat: PolyNTT = b_hat + a_vals
            with pytest.raises(TypeError):
                c_hat: PolyNTT = a_hat + b_vals[0]
            with pytest.raises(TypeError):
                c_hat: PolyNTT = a_vals[0] + b_hat
            with pytest.raises(TypeError):
                c_hat: PolyNTT = b_hat + a_vals[0]
            with pytest.raises(TypeError):
                c_hat: PolyNTT = a_hat + 0
            c_hat: PolyNTT = 0 + a_hat
            with pytest.raises(TypeError):
                c_hat: PolyNTT = a_hat + 1


def test_poly_ntt_sub():
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        root: int = find_prou(mod=q, deg=d)
        inv_root: int = pow(root, q - 2, q)
        for _ in range(SAMPLE_SIZE):
            a_vals: List[int] = [randrange(1, q) for _ in range(d)]
            a_hat: PolyNTT = PolyNTT(
                mod=q,
                vals=tuple(a_vals))
            b_vals: List[int] = [randrange(1, q) for _ in range(d)]
            b_hat: PolyNTT = PolyNTT(
                mod=q,
                vals=tuple(b_vals))
            c_hat: PolyNTT = a_hat - b_hat
            for i in range(d):
                assert (c_hat.vals[i] - (a_vals[i] - b_vals[i])) % q == 0
            with pytest.raises(TypeError):
                c_hat: PolyNTT = a_hat - b_vals
            with pytest.raises(TypeError):
                c_hat: PolyNTT = a_vals - b_hat
            with pytest.raises(TypeError):
                c_hat: PolyNTT = a_hat - b_vals
            with pytest.raises(TypeError):
                c_hat: PolyNTT = a_hat - b_vals[0]
            with pytest.raises(TypeError):
                c_hat: PolyNTT = a_vals[0] - b_hat
            with pytest.raises(TypeError):
                c_hat: PolyNTT = a_hat - 0
                assert c_hat == a_hat
            with pytest.raises(TypeError):
                c_hat: PolyNTT = 0 - a_hat
                assert c_hat == -a_hat
            with pytest.raises(TypeError):
                c_hat: PolyNTT = a_hat - 1


def test_poly_ntt_neg():
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        root: int = find_prou(mod=q, deg=d)
        inv_root: int = pow(root, q - 2, q)
        for _ in range(SAMPLE_SIZE):
            a_vals: List[int] = [randrange(1, q) for _ in range(d)]
            a_hat: PolyNTT = PolyNTT(
                mod=q,
                vals=tuple(a_vals))
            b_hat: PolyNTT = -a_hat
            for i in range(d):
                assert (b_hat.vals[i] + a_vals[i]) % q == 0


def test_poly_ntt_radd():
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        root: int = find_prou(mod=q, deg=d)
        inv_root: int = pow(root, q - 2, q)
        for _ in range(SAMPLE_SIZE):
            a_vals: List[int] = [randrange(1, q) for _ in range(d)]
            a_hat: PolyNTT = PolyNTT(
                mod=q,
                vals=tuple(a_vals))
            b_vals: List[int] = [randrange(1, q) for _ in range(d)]
            b_hat: PolyNTT = PolyNTT(
                mod=q,
                vals=tuple(b_vals))
            c_hat: PolyNTT = b_hat + a_hat
            for i in range(d):
                assert (c_hat.vals[i] - (a_vals[i] + b_vals[i])) % q == 0
            with pytest.raises(TypeError):
                c_hat: PolyNTT = b_hat + a_vals
            with pytest.raises(TypeError):
                c_hat: PolyNTT = b_vals + a_hat
            with pytest.raises(TypeError):
                c_hat: PolyNTT = b_hat + a_vals[0]
            with pytest.raises(TypeError):
                c_hat: PolyNTT = b_vals[0] + a_hat
            c_hat: PolyNTT = 0 + b_hat
            assert c_hat == b_hat
            with pytest.raises(TypeError):
                c_hat: PolyNTT = b_hat + q


def test_poly_ntt_mul():
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        root: int = find_prou(mod=q, deg=d)
        inv_root: int = pow(root, q - 2, q)
        for _ in range(SAMPLE_SIZE):
            a_vals: List[int] = [randrange(1, q) for _ in range(d)]
            a_hat: PolyNTT = PolyNTT(
                mod=q,
                vals=tuple(a_vals))
            b_vals: List[int] = [randrange(1, q) for _ in range(d)]
            b_hat: PolyNTT = PolyNTT(
                mod=q,
                vals=tuple(b_vals))
            c_hat: PolyNTT = a_hat * b_hat
            for i in range(d):
                assert (c_hat.vals[i] - (a_vals[i] * b_vals[i])) % q == 0
            with pytest.raises(TypeError):
                a_vals * b_hat
            with pytest.raises(TypeError):
                a_hat * b_vals
            with pytest.raises(TypeError):
                c_hat: PolyNTT = a_hat * 0
            assert c_hat == 0
            c_hat: int = 0 * b_hat
            assert c_hat == 0
            c_hat: PolyNTT = a_hat * 1
            assert c_hat == a_hat
            c_hat: PolyNTT = 1 * b_hat
            assert c_hat == b_hat


def test_poly_ntt_rmul():
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        root: int = find_prou(mod=q, deg=d)
        inv_root: int = pow(root, q - 2, q)
        for _ in range(SAMPLE_SIZE):
            a_vals: List[int] = [randrange(1, q) for _ in range(d)]
            a_hat: PolyNTT = PolyNTT(
                mod=q,
                vals=tuple(a_vals))
            b_vals: List[int] = [randrange(1, q) for _ in range(d)]
            b_hat: PolyNTT = PolyNTT(
                mod=q,
                vals=tuple(b_vals))
            c_hat: PolyNTT = b_hat * a_hat
            for i in range(d):
                assert (c_hat.vals[i] - (a_vals[i] * b_vals[i])) % q == 0
            with pytest.raises(TypeError):
                b_hat * a_vals
            with pytest.raises(TypeError):
                b_vals * a_hat
            with pytest.raises(TypeError):
                c_hat: PolyNTT = b_hat * 0
            assert c_hat == 0
            c_hat: int = 0 * a_hat
            assert c_hat == 0
            c_hat: PolyNTT = b_hat * 1
            assert c_hat == b_hat
            c_hat: PolyNTT = 1 * a_hat
            assert c_hat == a_hat


def test_transform_2d():
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        root: int = find_prou(mod=q, deg=d)
        inv_root: int = pow(root, q - 2, q)
        for _ in range(SAMPLE_SIZE):
            a_coefs: List[int] = [randrange(1, q) for _ in range(d)]
            a: Poly = Poly(
                mod=q,
                vals=tuple(a_coefs))
            b_coefs: List[int] = [randrange(1, q) for _ in range(d)]
            b: Poly = Poly(
                mod=q,
                vals=tuple(b_coefs))
            c_coefs: List[int] = [0 for _ in range(2 * d)]
            for i in range(len(a_coefs)):
                for j in range(len(b_coefs)):
                    c_coefs[i + j] = (c_coefs[i + j] + a_coefs[i] * b_coefs[j]) % q
            c_coefs = [(x - y) % q for x, y in zip(c_coefs[:d], c_coefs[d:])]
            c: Poly = Poly(
                mod=q,
                vals=tuple(c_coefs))

            a_ntt: PolyNTT = a.transform()
            b_ntt: PolyNTT = b.transform()
            a_ntt_times_b_ntt: PolyNTT = a_ntt * b_ntt
            inv_a_ntt_times_b_ntt: Poly = a_ntt_times_b_ntt.transform()
            assert inv_a_ntt_times_b_ntt == c


def test_comprehensive():
    for _ in range(SAMPLE_SIZE):
        for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
            root: int = find_prou(mod=q, deg=d)
            inv_root: int = pow(root, q - 2, q)
            # Random polynomial a
            a_coefs: List[int] = [randrange(1, q) for _ in range(d)]
            a: Poly = Poly(
                mod=q,
                vals=tuple(a_coefs))
            # Random polynomial b
            b_coefs: List[int] = [randrange(1, q) for _ in range(d)]
            b: Poly = Poly(
                mod=q,
                vals=tuple(b_coefs))

            # Transformed a and b
            a_hat: PolyNTT = a.transform()
            b_hat: PolyNTT = b.transform()

            # Product of transforms
            a_hat_times_b_hat: PolyNTT = a_hat * b_hat

            # Invert
            inv_a_hat_times_b_hat: Poly = a_hat_times_b_hat.transform()

            # Check
            assert inv_a_hat_times_b_hat == a * b
