import pytest
from api.ntt import ntt

TEST_DATA = [
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

@pytest.mark.parametrize("modulus, val, inv_flag, expected_result, observed_result", TEST_DATA)
def test_ntt_api(modulus, val, inv_flag, expected_result, observed_result):
    assert observed_result == expected_result
    z = ntt(val=val, mod=modulus, inv_flag=inv_flag)
    assert all((x-y) % modulus == 0 for x, y in zip(z, observed_result))