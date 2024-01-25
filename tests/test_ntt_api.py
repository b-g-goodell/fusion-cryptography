import pytest
from api.ntt_api import ntt_api

TEST_DATA = [
    (17, [1, 0, 0, 0, 0, 0, 0, 0], False, [1, 1, 1, 1, 1, 1, 1, 1], ntt_api(modulus=17, inv_flag=False, val=[1, 0, 0, 0, 0, 0, 0, 0])),
    (17, [0, 1, 0, 0, 0, 0, 0, 0], False, [3,-3,5,-5,-7,7,-6,6], ntt_api(modulus=17, inv_flag=False, val=[0, 1, 0, 0, 0, 0, 0, 0])),
    (17, [0, 0, 1, 0, 0, 0, 0, 0], False, [-8,-8,8,8,-2,-2,2,2], ntt_api(modulus=17, inv_flag=False, val=[0, 0, 1, 0, 0, 0, 0, 0])),
    (17, [0, 0, 0, 1, 0, 0, 0, 0], False, [-7,7,6,-6,-3,3,5,-5], ntt_api(modulus=17, inv_flag=False, val=[0, 0, 0, 1, 0, 0, 0, 0])),
    (17, [0, 0, 0, 0, 1, 0, 0, 0], False, [-4,-4,-4,-4,4,4,4,4], ntt_api(modulus=17, inv_flag=False, val=[0, 0, 0, 0, 1, 0, 0, 0])),
    (17, [0, 0, 0, 0, 0, 1, 0, 0], False, [5, -5, -3, 3, 6, -6, -7, 7], ntt_api(modulus=17, inv_flag=False, val=[0, 0, 0, 0, 0, 1, 0, 0])),
    (17, [0, 0, 0, 0, 0, 0, 1, 0], False, [-2, -2, 2, 2, -8, -8, 8, 8], ntt_api(modulus=17, inv_flag=False, val=[0, 0, 0, 0, 0, 0, 1, 0])),
    (17, [0, 0, 0, 0, 0, 0, 0, 1], False, [-6, 6, -7, 7, 5, -5, 3, -3], ntt_api(modulus=17, inv_flag=False, val=[0, 0, 0, 0, 0, 0, 0, 1])),
    (17, [1, 1, 1, 1, 1, 1, 1, 1], True, [1, 0, 0, 0, 0, 0, 0, 0], ntt_api(modulus=17, inv_flag=True, val=[1, 1, 1, 1, 1, 1, 1, 1])),
    (17, [3,-3,5,-5,-7,7,-6,6], True, [0, 1, 0, 0, 0, 0, 0, 0], ntt_api(modulus=17, inv_flag=True, val=[3,-3,5,-5,-7,7,-6,6])),
    (17, [-8,-8,8,8,-2,-2,2,2], True, [0, 0, 1, 0, 0, 0, 0, 0], ntt_api(modulus=17, inv_flag=True, val=[-8,-8,8,8,-2,-2,2,2])),
    (17, [-7,7,6,-6,-3,3,5,-5], True, [0, 0, 0, 1, 0, 0, 0, 0], ntt_api(modulus=17, inv_flag=True, val=[-7,7,6,-6,-3,3,5,-5])),
    (17, [-4,-4,-4,-4,4,4,4,4], True, [0, 0, 0, 0, 1, 0, 0, 0], ntt_api(modulus=17, inv_flag=True, val=[-4,-4,-4,-4,4,4,4,4])),
    (17, [5, -5, -3, 3, 6, -6, -7, 7], True, [0, 0, 0, 0, 0, 1, 0, 0], ntt_api(modulus=17, inv_flag=True, val=[5, -5, -3, 3, 6, -6, -7, 7])),
    (17, [-2, -2, 2, 2, -8, -8, 8, 8], True, [0, 0, 0, 0, 0, 0, 1, 0], ntt_api(modulus=17, inv_flag=True, val=[-2, -2, 2, 2, -8, -8, 8, 8])),
    (17, [-6, 6, -7, 7, 5, -5, 3, -3], True, [0, 0, 0, 0, 0, 0, 0, 1], ntt_api(modulus=17, inv_flag=True, val=[-6, 6, -7, 7, 5, -5, 3, -3])),
]

@pytest.mark.parametrize("modulus, val, inv_flag, expected_result, observed_result", TEST_DATA)
def test_ntt_api(modulus, val, inv_flag, expected_result, observed_result):
    assert observed_result == expected_result
    z = ntt_api(modulus=modulus, inv_flag=inv_flag, val=val)
    assert all((x-y) % modulus == 0 for x, y in zip(z, observed_result))