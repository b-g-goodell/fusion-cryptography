from secrets import randbits, randbelow
from typing import List
from api.polynomials import Polynomial as Poly
from api.matrices import Matrix
# from algebra.polynomials import _PolynomialCoefficientRepresentation, _PolynomialNTTRepresentation
# from algebra.matrices import _GeneralMatrix


def sample_polynomial_coefficient_representation(modulus: int, degree: int, norm_bound: int, weight_bound: int) -> Poly:
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
    return Poly(modulus=modulus, values=coefficients, representation='coefficient')


def sample_polynomial_ntt_representation(modulus: int, degree: int, root: int, inv_root: int, root_order: int) -> Poly:
    values: List[int] = [randbelow(modulus) - (modulus // 2) for _ in range(degree)]
    return Poly(modulus=modulus, values=values, representation='ntt')


def sample_coefficient_matrix(modulus: int, degree: int, root_order: int, root: int, inv_root: int, num_rows: int,
                              num_cols: int, norm_bound: int, weight_bound: int) -> Matrix:
    return Matrix(
        matrix=[
            [
                sample_polynomial_coefficient_representation(modulus=modulus, degree=degree, root=root,
                                                             inv_root=inv_root, root_order=root_order,
                                                             norm_bound=norm_bound, weight_bound=weight_bound)
                for j in range(num_cols)
            ]
            for i in range(num_rows)
        ]
    )


def sample_ntt_matrix(modulus: int, degree: int, root_order: int, root: int, inv_root: int, num_rows: int,
                      num_cols: int) -> Matrix:
    matrix: List[List[PolynomialNTTRepresentation]] = [[sample_polynomial_ntt_representation(modulus=modulus, degree=degree, root=root, inv_root=inv_root, root_order=root_order) for j in range(num_cols)] for i in range(num_rows)]
    return GeneralMatrix(matrix=matrix)

