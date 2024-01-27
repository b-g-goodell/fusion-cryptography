import pytest
from copy import deepcopy
from random import randrange
from typing import List, Tuple
# from algebra.ntt import find_prou, _ntt_poly_mult, _bit_reverse_copy
from api.ntt import ntt, find_prou, bit_reverse_copy, ntt_poly_mult
from algebra.polynomials import (
    _PolynomialCoefficientRepresentation as Poly,
    _PolynomialNTTRepresentation as PolyNTT,
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
            modulus=mod,
            degree=deg,
            root=root,
            inv_root=inv_root,
            values=a_coefs,
        )
        b: Poly = Poly(
            modulus=mod,
            degree=deg,
            root=root,
            inv_root=inv_root,
            values=b_coefs,
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
            modulus=mod,
            degree=deg,
            root=root,
            inv_root=inv_root,
            values=expected_c_coefs,
        )

        observed_c: Poly = a * b
        assert len(observed_c.values) == len(expected_c.values)
        assert all(
            (x - y) % mod == 0
            for x, y in zip(observed_c.values, expected_c.values)
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
        assert len(another_observed_c_coefs) == len(expected_c.values)
        assert all(
            (x - y) % mod == 0
            for x, y in zip(observed_c.values, another_observed_c_coefs))
        assert all(
            (x - y) % mod == 0
            for x, y in zip(expected_c.values, observed_c.values))
        another_observed_c: Poly = Poly(modulus=mod, degree=deg, root=root, inv_root=inv_root, values=another_observed_c_coefs)
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
                modulus=q,
                degree=d,
                root=root,
                inv_root=inv_root,
                values=a_coefs)

            b_coef: int = randrange(1, q)
            b_index: int = randrange(d)
            b_coefs: List[int] = [0 if i != b_index else b_coef for i in range(d)]
            b: Poly = Poly(
                modulus=q,
                degree=d,
                root=root,
                inv_root=inv_root,
                values=b_coefs)

            expected_c_coefs: List[int] = [0 for _ in range(d)]
            expected_c_coefs[(a_index + b_index) % d] = ((a_coef * b_coef) % q) * (
                1 - 2 * int(a_index + b_index >= d)
            )
            expected_c: Poly = Poly(
                modulus=q,
                degree=d,
                root=root,
                inv_root=inv_root,
                values=expected_c_coefs)
            observed_c: Poly = a * b
            try:
                assert expected_c == observed_c
            except AssertionError:
                expected_c: Poly = Poly(
                    modulus=q,
                    degree=d,
                    root=root,
                    inv_root=inv_root,
                    values=expected_c_coefs,
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
        Poly(modulus=1, degree=1, root=1, inv_root=1, values=1)
    with pytest.raises(TypeError):
        Poly(modulus=5, degree=2, root=1, inv_root=1, values=1)
    with pytest.raises(TypeError):
        Poly(
            modulus=1.0,
            degree=1,
            root=1,
            inv_root=1,
            values=["hello world"])
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        root: int = find_prou(mod=q, deg=d)
        inv_root: int = pow(root, q - 2, q)
        for _ in range(SAMPLE_SIZE):
            a_coefs: List[int] = [randrange(q) for _ in range(d)]
            a: Poly = Poly(
                modulus=q,
                degree=d,
                root=root,
                inv_root=inv_root,
                values=a_coefs)
            assert a.modulus == q
            assert a.degree == d
            assert a.root_order == 2 * d
            assert a.root == root
            assert a.inv_root == inv_root
            assert a.values == a_coefs


def test_poly_str():
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        root: int = find_prou(mod=q, deg=d)
        inv_root: int = pow(root, q - 2, q)
        for _ in range(SAMPLE_SIZE):
            a_coefs: List[int] = [randrange(q) for _ in range(d)]
            a: Poly = Poly(
                modulus=q,
                degree=d,
                root=root,
                inv_root=inv_root,
                values=a_coefs)
            assert (
                str(a)
                == f"PolynomialCoefficientRepresentation(modulus={q}, degree={d}, root={root}, inv_root={inv_root}, values={a_coefs})")


def test_poly_repr():
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        root: int = find_prou(mod=q, deg=d)
        inv_root: int = pow(root, q - 2, q)
        for _ in range(SAMPLE_SIZE):
            a_coefs: List[int] = [randrange(q) for _ in range(d)]
            a: Poly = Poly(
                modulus=q,
                degree=d,
                root=root,
                inv_root=inv_root,
                values=a_coefs)
            assert (
                repr(a)
                == f"PolynomialCoefficientRepresentation(modulus={q}, degree={d}, root={root}, inv_root={inv_root}, values={a_coefs})")


def test_poly_eq():
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        root: int = find_prou(mod=q, deg=d)
        inv_root: int = pow(root, q - 2, q)
        for _ in range(SAMPLE_SIZE):
            a_coefs: List[int] = [randrange(q) for _ in range(d)]
            a: Poly = Poly(
                modulus=q,
                degree=d,
                root=root,
                inv_root=inv_root,
                values=a_coefs)
            b: Poly = Poly(
                modulus=q,
                degree=d,
                root=root,
                inv_root=inv_root,
                values=a.values)
            assert a == a
            assert a == b
            assert len(a.values) == len(b.values)
            for i, next_coef in enumerate(a.values):
                a.values[i] = next_coef + randrange(2**10) * q
                b.values[i] = next_coef + randrange(2**10) * q
            assert a == b


def test_poly_add():
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        root: int = find_prou(mod=q, deg=d)
        inv_root: int = pow(root, q - 2, q)
        for _ in range(SAMPLE_SIZE):
            a_coefs: List[int] = [randrange(1, q) for _ in range(d)]
            b_coefs: List[int] = [randrange(1, q) for _ in range(d)]
            a: Poly = Poly(
                modulus=q,
                degree=d,
                root=root,
                inv_root=inv_root,
                values=a_coefs)
            b: Poly = Poly(
                modulus=q,
                degree=d,
                root=root,
                inv_root=inv_root,
                values=b_coefs)
            c: Poly = a + b
            assert c.modulus == q
            assert c.degree == d
            assert c.root_order == 2 * d
            assert c.root == root
            assert c.inv_root == inv_root
            for i in range(d):
                assert (c.values[i] - (a.values[i] + b.values[i])) % q == 0


def test_poly_sub():
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        root: int = find_prou(mod=q, deg=d)
        inv_root: int = pow(base=root, exp=q-2, mod=q)
        for _ in range(SAMPLE_SIZE):
            a_coefs: List[int] = [randrange(1, q) for _ in range(d)]
            b_coefs: List[int] = [randrange(1, q) for _ in range(d)]
            a: Poly = Poly(
                modulus=q,
                degree=d,
                root=root,
                inv_root=inv_root,
                values=a_coefs)
            b: Poly = Poly(
                modulus=q,
                degree=d,
                root=root,
                inv_root=inv_root,
                values=b_coefs)
            c: Poly = a - b
            assert c.modulus == q
            assert c.degree == d
            assert c.root_order == 2 * d
            assert c.root == root
            assert c.inv_root == inv_root
            for i in range(d):
                assert (c.values[i] - (a.values[i] - b.values[i])) % q == 0


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
                modulus=q,
                degree=d,
                root=root,
                inv_root=inv_root,
                values=a_coefs)
            b: Poly = Poly(
                modulus=q,
                degree=d,
                root=root,
                inv_root=inv_root,
                values=b_coefs)
            expected_c: Poly = Poly(
                modulus=q,
                degree=d,
                root=root,
                inv_root=inv_root,
                values=expected_c_coefs)
            observed_c: Poly = a * b
            assert expected_c == observed_c


def test_poly_norm():
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        root: int = find_prou(mod=q, deg=d)
        inv_root: int = pow(root, q - 2, q)
        for _ in range(SAMPLE_SIZE):
            a_coefs: List[int] = [randrange(1, q) for _ in range(d)]
            a: Poly = Poly(
                modulus=q,
                degree=d,
                root=root,
                inv_root=inv_root,
                values=a_coefs)
            expected_infinity_norm: int = max(abs(x) for x in a.values)
            observed_infinity_norm_and_weight: Tuple[int, int] = a.norm_weight
            assert expected_infinity_norm, expected_infinity_norm == observed_infinity_norm_and_weight


def test_poly_ntt_init():
    with pytest.raises(TypeError):
        PolyNTT(modulus=2, degree=1, root=1, inv_root=1, values=1)
    with pytest.raises(TypeError):
        PolyNTT(modulus=5, degree=2, root=-1, inv_root=-1, values=1)
    with pytest.raises(TypeError):
        PolyNTT(
            modulus=5,
            degree=2,
            root=1,
            inv_root=1,
            values=["hello world"])
    root = find_prou(mod=5, deg=2)
    inv_root = pow(root, 5 - 2, 5)
    with pytest.raises(TypeError):
        PolyNTT(
            modulus=5,
            degree=2,
            root=root,
            inv_root=inv_root,
            values=["hello world"])
    with pytest.raises(TypeError):
        PolyNTT(modulus=5, degree=2, root=root, inv_root=inv_root, values=[1])
    with pytest.raises(TypeError):
        PolyNTT(
            modulus=5,
            degree=2,
            root=root,
            inv_root=inv_root,
            values=[1, 2, 3])
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        root: int = find_prou(mod=q, deg=d)
        inv_root: int = pow(root, q - 2, q)
        for _ in range(SAMPLE_SIZE):
            a_coefs: List[int] = [randrange(1, q) for _ in range(d)]
            a: Poly = Poly(
                modulus=q,
                degree=d,
                root=root,
                inv_root=inv_root,
                values=a_coefs)
            b: Poly = Poly(
                modulus=q,
                degree=d,
                root=root,
                inv_root=inv_root,
                values=a.values)
            assert a == b
            assert a == a
            for i, next_coef in enumerate(a.values):
                a.values[i] = next_coef + randrange(2) * q
                b.values[i] = next_coef + randrange(2) * q
            assert a == b


def test_poly_ntt_str():
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        root: int = find_prou(mod=q, deg=d)
        inv_root: int = pow(root, q - 2, q)
        for _ in range(SAMPLE_SIZE):
            a_vals: List[int] = [randrange(1, q) for _ in range(d)]
            a_hat: PolyNTT = PolyNTT(
                modulus=q,
                degree=d,
                root=root,
                inv_root=inv_root,
                values=a_vals)
            assert (
                str(a_hat)
                == repr(a_hat)
                == f"PolynomialNTTRepresentation(modulus={q}, degree={d}, root={root}, inv_root={inv_root}, values={a_vals})")


def test_poly_ntt_eq():
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        root: int = find_prou(mod=q, deg=d)
        inv_root: int = pow(root, q - 2, q)
        for _ in range(SAMPLE_SIZE):
            a_vals: List[int] = [randrange(1, q) for _ in range(d)]
            a_hat: PolyNTT = PolyNTT(
                modulus=q,
                degree=d,
                root=root,
                inv_root=inv_root,
                values=a_vals)
            b_hat: PolyNTT = PolyNTT(
                modulus=q,
                degree=d,
                root=root,
                inv_root=inv_root,
                values=a_hat.values)
            assert a_hat == b_hat
            assert a_hat == a_hat
            for i, next_val in enumerate(a_hat.values):
                a_hat.values[i] = next_val + randrange(2) * q
                b_hat.values[i] = next_val + randrange(2) * q
            assert a_hat == b_hat


def test_poly_ntt_add():
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        root: int = find_prou(mod=q, deg=d)
        inv_root: int = pow(root, q - 2, q)
        for _ in range(SAMPLE_SIZE):
            a_vals: List[int] = [randrange(1, q) for _ in range(d)]
            a_hat: PolyNTT = PolyNTT(
                modulus=q,
                degree=d,
                root=root,
                inv_root=inv_root,
                values=a_vals)
            b_vals: List[int] = [randrange(1, q) for _ in range(d)]
            b_hat: PolyNTT = PolyNTT(
                modulus=q,
                degree=d,
                root=root,
                inv_root=inv_root,
                values=b_vals)
            c_hat: PolyNTT = a_hat + b_hat
            assert all([
                (z - (x + y)) % q == 0 for x, y, z in zip(a_vals, b_vals, c_hat.values)
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
                modulus=q,
                degree=d,
                root=root,
                inv_root=inv_root,
                values=a_vals)
            b_vals: List[int] = [randrange(1, q) for _ in range(d)]
            b_hat: PolyNTT = PolyNTT(
                modulus=q,
                degree=d,
                root=root,
                inv_root=inv_root,
                values=b_vals)
            c_hat: PolyNTT = a_hat - b_hat
            for i in range(d):
                assert (c_hat.values[i] - (a_vals[i] - b_vals[i])) % q == 0
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
            c_hat: PolyNTT = a_hat - 0
            with pytest.raises(TypeError):
                c_hat: PolyNTT = 0 - a_hat
            with pytest.raises(TypeError):
                c_hat: PolyNTT = a_hat - 1


def test_poly_ntt_neg():
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        root: int = find_prou(mod=q, deg=d)
        inv_root: int = pow(root, q - 2, q)
        for _ in range(SAMPLE_SIZE):
            a_vals: List[int] = [randrange(1, q) for _ in range(d)]
            a_hat: PolyNTT = PolyNTT(
                modulus=q,
                degree=d,
                root=root,
                inv_root=inv_root,
                values=a_vals)
            b_hat: PolyNTT = -a_hat
            for i in range(d):
                assert (b_hat.values[i] + a_vals[i]) % q == 0


def test_poly_ntt_radd():
    for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
        root: int = find_prou(mod=q, deg=d)
        inv_root: int = pow(root, q - 2, q)
        for _ in range(SAMPLE_SIZE):
            a_vals: List[int] = [randrange(1, q) for _ in range(d)]
            a_hat: PolyNTT = PolyNTT(
                modulus=q,
                degree=d,
                root=root,
                inv_root=inv_root,
                values=a_vals)
            b_vals: List[int] = [randrange(1, q) for _ in range(d)]
            b_hat: PolyNTT = PolyNTT(
                modulus=q,
                degree=d,
                root=root,
                inv_root=inv_root,
                values=b_vals)
            c_hat: PolyNTT = b_hat + a_hat
            for i in range(d):
                assert (c_hat.values[i] - (a_vals[i] + b_vals[i])) % q == 0
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
                modulus=q,
                degree=d,
                root=root,
                inv_root=inv_root,
                values=a_vals)
            b_vals: List[int] = [randrange(1, q) for _ in range(d)]
            b_hat: PolyNTT = PolyNTT(
                modulus=q,
                degree=d,
                root=root,
                inv_root=inv_root,
                values=b_vals)
            c_hat: PolyNTT = a_hat * b_hat
            for i in range(d):
                assert (c_hat.values[i] - (a_vals[i] * b_vals[i])) % q == 0
            with pytest.raises(TypeError):
                a_hat * b_vals
            with pytest.raises(TypeError):
                a_vals * b_hat

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
                modulus=q,
                degree=d,
                root=root,
                inv_root=inv_root,
                values=a_vals)
            b_vals: List[int] = [randrange(1, q) for _ in range(d)]
            b_hat: PolyNTT = PolyNTT(
                modulus=q,
                degree=d,
                root=root,
                inv_root=inv_root,
                values=b_vals)
            c_hat: PolyNTT = b_hat * a_hat
            for i in range(d):
                assert (c_hat.values[i] - (a_vals[i] * b_vals[i])) % q == 0
            with pytest.raises(TypeError):
                b_hat * a_vals
            with pytest.raises(TypeError):
                b_vals * a_hat

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
                modulus=q,
                degree=d,
                root=root,
                inv_root=inv_root,
                values=a_coefs)
            a_hat: PolyNTT = a.transform()

            b_coefs: List[int] = [randrange(1, q) for _ in range(d)]
            b: Poly = Poly(
                modulus=q,
                degree=d,
                root=root,
                inv_root=inv_root,
                values=b_coefs)
            b_hat: PolyNTT = b.transform()

            c_coefs: List[int] = [0 for _ in range(2 * d)]
            for i in range(len(a_coefs)):
                for j in range(len(b_coefs)):
                    c_coefs[i + j] = (c_coefs[i + j] + a_coefs[i] * b_coefs[j]) % q
            c_coefs = [(x - y) % q for x, y in zip(c_coefs[:d], c_coefs[d:])]
            c: Poly = Poly(
                modulus=q,
                degree=d,
                root=root,
                inv_root=inv_root,
                values=c_coefs)

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
                modulus=q,
                degree=d,
                root=root,
                inv_root=inv_root,
                values=a_coefs)
            # Random polynomial b
            b_coefs: List[int] = [randrange(1, q) for _ in range(d)]
            b: Poly = Poly(
                modulus=q,
                degree=d,
                root=root,
                inv_root=inv_root,
                values=b_coefs)

            # Transformed a and b
            a_hat: PolyNTT = a.transform()
            b_hat: PolyNTT = b.transform()

            # Product of transforms
            a_hat_times_b_hat: PolyNTT = a_hat * b_hat

            # Invert
            inv_a_hat_times_b_hat: Poly = a_hat_times_b_hat.transform()

            # Check
            assert inv_a_hat_times_b_hat == a * b

#
# def test_sample_polynomial_coefficient_representation():
#     modulus: int = 65537
#     degree: int = 1024
#     root_order: int = 2 * degree
#     root: int = find_prou(mod=modulus, deg=degree)
#     inv_root: int = pow(root, modulus - 2, modulus)
#     norm_bound: int = 1000
#     weight_bound: int = 100
#     seed: int = 123456789
#     f: Poly = sample_polynomial_coefficient_representation(modulus=modulus, degree=degree, root=root, inv_root=inv_root,
#                                                            root_order=root_order, norm_bound=norm_bound,
#                                                            weight_bound=weight_bound)
#     assert isinstance(f, Poly)
#     assert f.modulus == modulus
#     assert f.degree == degree
#     assert f.root == root
#     assert f.inv_root == inv_root
#     assert f.root_order == root_order
#     assert len(f.values) == degree
#     norm, weight = f.norm_weight
#     assert norm <= norm_bound
#     assert weight <= weight_bound
