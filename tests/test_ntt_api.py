import pytest
from api.ntt_api import ntt_api

TEST_DATA = [
    (17, [1, 0, 0, 0, 0, 0, 0, 0], [1, 1, 1, 1, 1, 1, 1, 1],     ntt_api(modulus=17, inv_flag=False, val=[1, 0, 0, 0, 0, 0, 0, 0]),     ntt_api(modulus=17, inv_flag=True,  val=[1, 1, 1, 1, 1, 1, 1, 1])),
    (17, [0, 1, 0, 0, 0, 0, 0, 0], [-8, 8, 8, -8, -2, 2, 2, -2], ntt_api(modulus=17, inv_flag=False, val=[0, 1, 0, 0, 0, 0, 0, 0]),     ntt_api(modulus=17, inv_flag=True,  val=[-8, 8, 8, -8, -2, 2, 2, -2])),
    (17, [0, 0, 1, 0, 0, 0, 0, 0], [-4, -4, 4, 4, 4, 4, -4, -4], ntt_api(modulus=17, inv_flag=False, val=[0, 0, 1, 0, 0, 0, 0, 0]),     ntt_api(modulus=17, inv_flag=True,  val=[-4, -4, 4, 4, 4, 4, -4, -4])),
    (17, [0, 0, 0, 1, 0, 0, 0, 0], [-2, 2, -2, 2, -8, 8, -8, 8], ntt_api(modulus=17, inv_flag=False, val=[0, 0, 0, 1, 0, 0, 0, 0]),     ntt_api(modulus=17, inv_flag=True,  val=[-2, 2, -2, 2, -8, 8, -8, 8])),
    (17, [0, 0, 0, 0, 1, 0, 0, 0], [-1, -1, -1, -1, 1, 1, 1, 1], ntt_api(modulus=17, inv_flag=False, val=[0, 0, 0, 0, 1, 0, 0, 0]),     ntt_api(modulus=17, inv_flag=True,  val=[-1, -1, -1, -1, 1, 1, 1, 1])),
    (17, [0, 0, 0, 0, 0, 1, 0, 0], [8, -8, -8, 8, -2, 2, 2, -2], ntt_api(modulus=17, inv_flag=False, val=[0, 0, 0, 0, 0, 1, 0, 0]),     ntt_api(modulus=17, inv_flag=True,  val=[8, -8, -8, 8, -2, 2, 2, -2])),
    (17, [0, 0, 0, 0, 0, 0, 1, 0], [4, 4, -4, -4, 4, 4, -4, -4], ntt_api(modulus=17, inv_flag=False, val=[0, 0, 0, 0, 0, 0, 1, 0]),     ntt_api(modulus=17, inv_flag=True,  val=[4, 4, -4, -4, 4, 4, -4, -4])),
    (17, [0, 0, 0, 0, 0, 0, 0, 1], [2, -2, 2, -2, -8, 8, -8, 8], ntt_api(modulus=17, inv_flag=False, val=[0, 0, 0, 0, 0, 0, 0, 1]),     ntt_api(modulus=17, inv_flag=True,  val=[2, -2, 2, -2, -8, 8, -8, 8])),
]

@pytest.mark.parametrize("modulus, a, expected_a_hat, observed_a_hat, observed_a_hat_hat", TEST_DATA)
def test_ntt_api(modulus, a, expected_a_hat, observed_a_hat, observed_a_hat_hat):
    assert all((x-y) % modulus == 0 for x, y in zip(expected_a_hat, observed_a_hat))
    assert all((x-y) % modulus == 0 for x, y in zip(a, observed_a_hat_hat))