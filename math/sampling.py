from secrets import randbits, randbelow
from typing import Dict, List, Union, Optional
from math.polynomials import PolynomialCoefficientRepresentation, PolynomialNTTRepresentation
def sample_polynomial_coefficient_representation(modulus: int, degree: int, root: int, inv_root: int, root_order: int, norm_bound: int, weight_bound: int) -> PolynomialCoefficientRepresentation:
    # Exactly weight non-zero coefficients
    num_coefs_to_gen: int = max(0, min(degree, weight_bound))
    bound: int = max(0, min(modulus // 2, norm_bound))
    coefficients: List[int] = [
        (1 + randbelow(bound)) * (1 - 2 * randbits(1)) for _ in range(num_coefs_to_gen)
    ]
    coefficients += [0 for _ in range(degree - len(coefficients))]
    if num_coefs_to_gen < degree:
        # fisher-yates shuffle
        for i in range(degree - 1, 0, -1):
            j = randbelow(i + 1)
            coefficients[i], coefficients[j] = coefficients[j], coefficients[i]
    return PolynomialCoefficientRepresentation(
        modulus=modulus,
        degree=degree,
        root=root,
        inv_root=inv_root,
        root_order=root_order,
        coefficients=coefficients,
    )


def sample_polynomial_ntt_representation(modulus: int, degree: int, root: int, inv_root: int, root_order: int) -> PolynomialNTTRepresentation:
    values: List[int] = [randbelow(modulus) - (modulus // 2) for _ in range(degree)]
    return PolynomialNTTRepresentation(
        modulus=modulus,
        degree=degree,
        root=root,
        inv_root=inv_root,
        root_order=root_order,
        values=values,
    )
