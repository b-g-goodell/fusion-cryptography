# import pytest
# from typing import List, Dict
# from algebra.matrices import GeneralMatrix
# from algebra.sampling import sample_coefficient_matrix, sample_ntt_matrix
# from fusion.fusion import PREFIX_PARAMETERS, Params, fusion_setup  # for dummy testing values
#
#
# def test_prefix_parameters():
#     assert 128 in PREFIX_PARAMETERS
#     assert 256 in PREFIX_PARAMETERS
#
#
# def test_sample_coefficient_matrix():
#     # Setup
#     params: Params = fusion_setup(secpar=128)
#
#     x: GeneralMatrix = sample_coefficient_matrix(modulus=params.modulus, degree=params.degree,
#                                                  root_order=params.root_order, root=params.root,
#                                                  inv_root=params.inv_root, num_rows=1, num_cols=1, norm_bound=1,
#                                                  weight_bound=1)
#     assert len(x.matrix) == 1
#     assert all([len(row) == 1 for row in x.matrix])
#     assert all(hasattr(z, "norm") for y in x.matrix for z in y)
#     assert all(hasattr(z, "weight") for y in x.matrix for z in y)
#     assert 0 <= x.norm(p="infty") <= 1
#     assert 0 <= x.weight() <= 1
#
#     x: GeneralMatrix = sample_coefficient_matrix(modulus=params.modulus, degree=params.degree,
#                                                  root_order=params.root_order, root=params.root,
#                                                  inv_root=params.inv_root, num_rows=2, num_cols=3, norm_bound=17,
#                                                  weight_bound=16)
#     assert len(x.matrix) == 2
#     assert all([len(row) == 3 for row in x.matrix])
#     assert all(hasattr(z, "norm") for y in x.matrix for z in y)
#     assert all(hasattr(z, "weight") for y in x.matrix for z in y)
#     assert 0 <= x.norm(p="infty") <= 17
#     assert 0 <= x.weight() <= 16
#
#
# # @pytest.mark.skip
# def test_sample_ntt_matrix():
#     # Setup
#     params: Params = fusion_setup(secpar=128)
#
#     x: GeneralMatrix = sample_ntt_matrix(modulus=params.modulus, degree=params.degree, root_order=params.root_order,
#                                          root=params.root, inv_root=params.inv_root, num_rows=1, num_cols=1)
#     assert len(x.matrix) == 1
#     assert all([len(row) == 1 for row in x.matrix])
#
#     x: GeneralMatrix = sample_ntt_matrix(modulus=params.modulus, degree=params.degree, root_order=params.root_order,
#                                          root=params.root, inv_root=params.inv_root, num_rows=2, num_cols=3)
#     assert len(x.matrix) == 2
#     assert all([len(row) == 3 for row in x.matrix])
