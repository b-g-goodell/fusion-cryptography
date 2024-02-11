from secrets import randbits, randbelow
from typing import List
from api.polynomials import Polynomial as Poly
from api.matrices import Matrix


def sample_poly_by_coefs(mod: int, deg: int, norm: int, wght: int) -> Poly:
    # Exactly weight non-zero coefficients
    num_coefs_to_gen: int = max(0, min(deg, wght))
    bound: int = max(0, min(mod // 2, norm))
    coefficients: List[int] = [
        (1 + randbelow(bound)) * (1 - 2 * randbits(1)) for _ in range(num_coefs_to_gen)
    ]
    coefficients += [0 for _ in range(deg - len(coefficients))]
    if num_coefs_to_gen < deg:
        # fisher-yates shuffle
        for i in range(deg - 1, 0, -1):
            j = randbelow(i + 1)
            coefficients[i], coefficients[j] = coefficients[j], coefficients[i]
    return Poly(mod=mod, vals=coefficients, rep_flag='coefficient')


def sample_poly_by_ntt(mod: int, deg: int) -> Poly:
    values: List[int] = [randbelow(mod) - (mod // 2) for _ in range(deg)]
    return Poly(mod=mod, vals=values, rep_flag='ntt')


def sample_matrix_by_coefs(mod: int, deg: int, num_rows: int, num_cols: int, norm: int, wght: int) -> Matrix:
    return Matrix(matrix=[
        [sample_poly_by_coefs(mod=mod, deg=deg, norm=norm, wght=wght) for j in range(num_cols)]
        for i in range(num_rows)])


def sample_matrix_by_ntt(modulus: int, degree: int, num_rows: int, num_cols: int) -> Matrix:
    matrix: List[List[Poly]] = [[sample_poly_by_ntt(mod=modulus, deg=degree)
                                 for j in range(num_cols)] for i in range(num_rows)]
    return Matrix(matrix=matrix)

