import pytest
from copy import deepcopy
from random import randrange
from typing import List, Tuple
from algebra.polynomials import (
    PolynomialCoefficientRepresentation as Poly,
    PolynomialNTTRepresentation as PolyNTT,
    transform,
)
from algebra.ntt import cent, bit_reverse_copy, ntt_poly_mult, cooley_tukey_ntt, gentleman_sande_intt, find_primitive_root
from test_ntt import SAMPLE_SIZE, PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE


@pytest.fixture
@pytest.mark.parametrize("d, q", PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE)
def one(d, q):
    return Poly(modulus=q, representation=[1] + [0] * (d-1))


@pytest.fixture
def ntt_one(one):
    return transform(one)


@pytest.mark.parametrize("d, q", PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE)
def test_polynomial_ntt_consistency(d: int, q: int):
    root = find_primitive_root(modulus=q, root_order=2*d)
    root_powers = [pow(base=root, exp=i, mod=q) for i in range(2*d)]
    inv_root = pow(base=root, exp=q-2, mod=q)
    inv_root_powers = [pow(base=inv_root, exp=i, mod=q) for i in range(2*d)]
    brv_root_powers = bit_reverse_copy(val=root_powers)
    brv_inv_root_powers = bit_reverse_copy(val=inv_root_powers)
    halfmod = q//2
    logmod = q.bit_length()

    # Create two random polynomials f and g
    coefs_f = [randrange(q) for _ in range(d)]
    coefs_g = [randrange(q) for _ in range(d)]
    expected_fg = [0 for _ in range(2*d)]
    for i, x in enumerate(coefs_f):
        for j, y in enumerate(coefs_g):
            expected_fg[i+j] += x*y
    # expected_fg = [cent(val=x, modulus=q, halfmod=halfmod, logmod=logmod) for x in expected_fg]
    expected_fg = [cent(val=x-y, modulus=q, halfmod=halfmod, logmod=logmod) for x, y in zip(expected_fg[:d], expected_fg[d:])]

    # Using polynomials.py for standard domain multiplication
    f = Poly(modulus=q, representation=coefs_f)
    g = Poly(modulus=q, representation=coefs_g)
    fg_standard = f * g

    coefs_fg_standard = fg_standard.representation
    assert all((x-y) % q == 0 for x, y in zip(coefs_fg_standard, expected_fg))

    # Using ntt.py for NTT domain multiplication
    f_ntt = cooley_tukey_ntt(val=deepcopy(coefs_f), modulus=q, root_order=2*d, bit_rev_root_powers=brv_root_powers, halfmod=halfmod, logmod=logmod)
    g_ntt = cooley_tukey_ntt(val=deepcopy(coefs_g), modulus=q, root_order=2*d, bit_rev_root_powers=brv_root_powers, halfmod=halfmod, logmod=logmod)
    fg_ntt = [cent(val=x*y, modulus=q, halfmod=halfmod, logmod=logmod) for x, y in zip(f_ntt, g_ntt)]
    fg_standard_ntt = cooley_tukey_ntt(val=deepcopy(coefs_fg_standard), modulus=q, root_order=2*d, bit_rev_root_powers=brv_root_powers, halfmod=halfmod, logmod=logmod)

    fg_reconstructed = gentleman_sande_intt(val=deepcopy(fg_ntt), modulus=q, root_order=2*d, bit_rev_inv_root_powers=brv_inv_root_powers, halfmod=halfmod, logmod=logmod)
    fg_standard_reconstructed = gentleman_sande_intt(val=deepcopy(fg_standard_ntt), modulus=q, root_order=2*d, bit_rev_inv_root_powers=brv_inv_root_powers, halfmod=halfmod, logmod=logmod)

    # Compare results
    assert all((x-y) % q == 0 for x, y in zip(fg_standard_ntt, fg_ntt))


PRODUCT_TEST_DATA: List[Tuple] = []
for d, q in PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE:
    for i in range(0, d, d//4):
        for j in range(0, d, d//4):
            coefs_f: List[int] = [randrange(q) if k == i else 0 for k in range(d)]
            coefs_g: List[int] = [randrange(q) if k == j else 0 for k in range(d)]
            coefs_fg: List[int] = [(sum(coefs_f)*sum(coefs_g)*(-1 if i + j >= d else 1)) % q if (k - (i+j)) % q == 0 else 0 for k in range(d)]
            f: Poly = Poly(modulus=q, representation=coefs_f)
            g: Poly = Poly(modulus=q, representation=coefs_g)
            fg: Poly = Poly(modulus=q, representation=coefs_fg)
            f_hat: PolyNTT = transform(f)
            g_hat: PolyNTT = transform(g)
            fg_hat: PolyNTT = transform(fg)
            f_hat_times_g_hat: PolyNTT = f_hat*g_hat
            f_hat_hat: Poly = transform(f_hat)
            g_hat_hat: Poly = transform(g_hat)
            f_hat_times_g_hat_hat: Poly = transform(f_hat_times_g_hat)
            fg_hat_hat: Poly = transform(fg_hat)
            PRODUCT_TEST_DATA += [tuple([
                d, q, i, j, sum(coefs_f), sum(coefs_g), sum(coefs_fg), coefs_f, coefs_g, coefs_fg, f, g, fg, f_hat,
                g_hat, fg_hat, f_hat_times_g_hat, f_hat_hat, g_hat_hat, f_hat_times_g_hat_hat, fg_hat_hat
            ])]


@pytest.mark.parametrize("d, q", PAIRS_OF_D_AND_Q_FORCING_ROU_EXISTENCE)
def test_ntt_polynomial_multiplication_failing_case(d: int, q: int):
    i, j = randrange(d), randrange(d)
    coefs_f = [1+randrange(q-1) if k == i else 0 for k in range(d)]
    coefs_g = [1+randrange(q-1) if k == j else 0 for k in range(d)]

    f = Poly(modulus=q, representation=coefs_f)
    g = Poly(modulus=q, representation=coefs_g)

    fg = f * g

    f_hat = transform(f)
    g_hat = transform(g)
    fg_hat = transform(fg)

    f_hat_times_g_hat = f_hat * g_hat

    fg_reconstructed = transform(f_hat_times_g_hat)

    assert fg_hat == f_hat_times_g_hat
    assert fg == fg_reconstructed


@pytest.mark.parametrize("d,q,i,j,const1,const2,const1const2,coefs_f,coefs_g,coefs_fg,f,g,fg,f_hat,g_hat,fg_hat,"
                         + "f_hat_times_g_hat,f_hat_hat,g_hat_hat,f_hat_times_g_hat_hat,fg_hat_hat", PRODUCT_TEST_DATA)
def test_polynomials(d: int, q: int, i: int, j: int, const1: int, const2: int, const1const2: int, coefs_f: List[int],
                     coefs_g: List[int], coefs_fg: List[int], f: Poly, g: Poly, fg: Poly, f_hat: PolyNTT,
                     g_hat: PolyNTT, fg_hat: PolyNTT, f_hat_times_g_hat: PolyNTT, f_hat_hat: Poly, g_hat_hat: Poly,
                     f_hat_times_g_hat_hat: Poly, fg_hat_hat: Poly):
    assert const1 == sum(coefs_f)
    assert const2 == sum(coefs_g)
    assert (abs(const1*const2) - abs(const1const2)) % q == 0
    assert (abs(const1const2) - abs(sum(coefs_fg))) % q == 0
    if i + j >= d:
        assert (sum(coefs_fg) + const1*const2) % q == 0
    else:
        assert (sum(coefs_fg) - const1*const2) % q == 0
    assert f*g == fg
    assert f_hat*g_hat == f_hat_times_g_hat
    assert fg_hat == transform(fg)
    assert fg_hat == f_hat_times_g_hat
    assert fg == f_hat_times_g_hat_hat
