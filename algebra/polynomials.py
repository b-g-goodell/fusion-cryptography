from copy import deepcopy
from typing import Dict, List, Union
from algebra.ntt import (
    _cent,
    _cooley_tukey_ntt,
    _gentleman_sande_intt,
    _bit_reverse_copy,
    _has_prou,
    _is_prou,
    _is_root_inverse,
)
from algebra.errors import *

cached_halfmods: Dict[int, int] = {}
cached_logmods: Dict[int, int] = {}


class PolynomialRepresentation(object):
    modulus: int
    degree: int
    root: int
    inv_root: int
    root_order: int

    def __init__(
        self, modulus: int, degree: int, root: int, inv_root: int, root_order: int
    ):
        if not isinstance(modulus, int) or not isinstance(degree, int) or not isinstance(root, int) or not isinstance(inv_root, int) or not isinstance(root_order, int):
            raise TypeError(MUST_BE_INT_ERR)
        elif not _has_prou(mod=modulus, deg=degree):
            raise ValueError(MUST_HAVE_PROU_ERR)
        elif not _is_prou(root=root, mod=modulus, deg=degree):
            raise ValueError(MUST_BE_CORRECT_ROOT_ERR)
        elif not _is_root_inverse:
            raise ValueError(MUST_BE_CORRECT_INVERSE_ROOT_ERR)
        self.modulus = modulus
        self.degree = degree
        self.root = root
        self.inv_root = inv_root
        self.root_order = root_order

    @property
    def halfmod(self) -> int:
        if self.modulus not in cached_halfmods:
            cached_halfmods[self.modulus] = self.modulus // 2
        return cached_halfmods[self.modulus]

    @property
    def logmod(self) -> int:
        if self.modulus not in cached_logmods:
            cached_logmods[self.modulus] = self.modulus.bit_length()
        return cached_logmods[self.modulus]


class PolynomialCoefficientRepresentation(PolynomialRepresentation):
    coefficients: List[int]

    def __init__(
        self,
        modulus: int,
        degree: int,
        root: int,
        inv_root: int,
        root_order: int,
        coefficients: List[int],
    ):
        super().__init__(
            modulus=modulus,
            degree=degree,
            root=root,
            inv_root=inv_root,
            root_order=root_order,
        )
        if not isinstance(coefficients, list):
            raise TypeError(MUST_BE_LIST_ERR)
        elif not all(isinstance(x, int) for x in coefficients):
            raise TypeError(MUST_BE_INT_ERR)
        elif len(coefficients) != degree:
            raise ValueError(DEGREE_MISMATCH_ERR)
        self.coefficients = coefficients

    def __str__(self):
        return f"PolynomialCoefficientRepresentation(modulus={self.modulus}, degree={self.degree}, root={self.root}, inv_root={self.inv_root}, root_order={self.root_order}, coefficients={self.coefficients})"

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        if not isinstance(other, PolynomialCoefficientRepresentation):
            return False
        elif self.modulus != other.modulus:
            return False
        elif self.degree != other.degree:
            return False
        elif self.root != other.root:
            return False
        elif self.root_order != other.root_order:
            return False
        return all(
            (x - y) % self.modulus == 0
            for x, y in zip(self.coefficients, other.coefficients)
        )

    def __add__(self, other):
        if other == 0:
            return self
        elif not isinstance(other, PolynomialCoefficientRepresentation):
            raise NotImplementedError(TYPE_MISMATCH_ERR)
        elif self.modulus != other.modulus:
            raise NotImplementedError(MODULUS_MISMATCH_ERR)
        elif self.degree != other.degree:
            raise NotImplementedError(DEGREE_MISMATCH_ERR)
        elif self.root_order != other.root_order or self.root_order != 2*self.degree:
            raise NotImplementedError(ROOT_ORDER_MISMATCH_ERR)
        elif self.root != other.root:
            raise NotImplementedError(ROOT_MISMATCH_ERR)
        return PolynomialCoefficientRepresentation(
            modulus=self.modulus,
            degree=self.degree,
            root=self.root,
            inv_root=self.inv_root,
            root_order=self.root_order,
            coefficients=[
                _cent(val=x + y, mod=self.modulus, halfmod=self.halfmod, logmod=self.logmod)
                for x, y in zip(self.coefficients, other.coefficients)
            ],
        )

    def __radd__(self, other):
        if other == 0:
            return self
        return self + other

    def __neg__(self):
        return PolynomialCoefficientRepresentation(
            modulus=self.modulus,
            degree=self.degree,
            root=self.root,
            inv_root=self.inv_root,
            root_order=self.root_order,
            coefficients=[-(x % self.modulus) for x in self.coefficients],
        )

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return other + (-self)

    def __mul__(self, other):
        if other == 0:
            return 0
        elif other == 1:
            return self
        elif not isinstance(other, PolynomialCoefficientRepresentation):
            raise NotImplementedError(TYPE_MISMATCH_ERR)
        elif self.modulus != other.modulus:
            raise NotImplementedError(MODULUS_MISMATCH_ERR)
        elif self.degree != other.degree:
            raise NotImplementedError(DEGREE_MISMATCH_ERR)
        elif self.root_order != other.root_order or self.root_order != 2*self.degree:
            raise NotImplementedError(ROOT_ORDER_MISMATCH_ERR)
        elif self.root != other.root:
            raise NotImplementedError(ROOT_MISMATCH_ERR)
        c: List[int] = [0 for _ in range(2 * self.degree)]
        for i, x in enumerate(self.coefficients):
            for j, y in enumerate(other.coefficients):
                c[i + j] += x * y
        c = [_cent(val=x - y, mod=self.modulus, halfmod=self.halfmod, logmod=self.logmod)
             for x, y in zip(c[:self.degree], c[self.degree:])]
        return PolynomialCoefficientRepresentation(
            modulus=self.modulus,
            degree=self.degree,
            root=self.root,
            inv_root=self.inv_root,
            root_order=self.root_order,
            coefficients=c,
        )

    def __rmul__(self, other):
        return self.__mul__(other=other)

    def norm(self, p: Union[int, str]) -> int:
        if p != "infty":
            raise NotImplementedError(NORM_P_NOT_IMPLEMENTED_ERR)
        return max(abs(x) for x in self.coefficients)

    def weight(self) -> int:
        return sum(1 if x % self.modulus != 0 else 0 for x in self.coefficients)


class PolynomialNTTRepresentation(PolynomialRepresentation):
    values: List[int]

    def __init__(
        self,
        modulus: int,
        degree: int,
        root: int,
        inv_root: int,
        root_order: int,
        values: List[int],
    ):
        super().__init__(
            modulus=modulus,
            degree=degree,
            root=root,
            inv_root=inv_root,
            root_order=root_order,
        )
        if not isinstance(values, list):
            raise TypeError(MUST_BE_LIST_ERR)
        elif not all(isinstance(x, int) for x in values):
            raise TypeError(MUST_BE_INT_ERR)
        elif len(values) != degree:
            raise ValueError(DEGREE_MISMATCH_ERR)
        self.values = values

    def __str__(self):
        return f"PolynomialNTTRepresentation(modulus={self.modulus}, degree={self.degree}, root={self.root}, inv_root={self.inv_root}, root_order={self.root_order}, values={self.values})"

    def __repr__(self):
        return self.__str__()

    def __eq__(self, other):
        if isinstance(other, int) and other == 0:
            return all(x % self.modulus == 0 for x in self.values)
        elif not isinstance(other, PolynomialNTTRepresentation):
            return False
        elif self.modulus != other.modulus:
            return False
        elif self.degree != other.degree:
            return False
        elif self.root_order != other.root_order:
            return False
        elif self.root != other.root or self.inv_root != other.inv_root:
            return False
        elif len(self.values) != len(other.values):
            return False
        return all(
            (x - y) % self.modulus == 0 for x, y in zip(self.values, other.values)
        )

    def __add__(self, other):
        if other == 0:
            return self
        elif not isinstance(other, PolynomialNTTRepresentation):
            raise NotImplementedError(TYPE_MISMATCH_ERR)
        elif self.modulus != other.modulus:
            raise NotImplementedError(MODULUS_MISMATCH_ERR)
        elif self.degree != other.degree or len(self.values) != self.degree or other.degree != len(other.values):
            raise NotImplementedError(DEGREE_MISMATCH_ERR)
        elif self.root_order != other.root_order or self.root_order != 2*self.degree or other.root_order != 2*other.degree:
            raise NotImplementedError(ROOT_ORDER_MISMATCH_ERR)
        elif self.root != other.root:
            raise NotImplementedError(ROOT_MISMATCH_ERR)
        return PolynomialNTTRepresentation(
            modulus=self.modulus,
            degree=self.degree,
            root=self.root,
            inv_root=self.inv_root,
            root_order=self.root_order,
            values=[
                _cent(val=x + y, mod=self.modulus, halfmod=self.halfmod, logmod=self.logmod)
                for x, y in zip(self.values, other.values)
            ],
        )

    def __radd__(self, other):
        if other == 0:
            return self
        return self + other

    def __neg__(self):
        return PolynomialNTTRepresentation(
            modulus=self.modulus,
            degree=self.degree,
            root=self.root,
            inv_root=self.inv_root,
            root_order=self.root_order,
            values=[-(x % self.modulus) for x in self.values],
        )

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return other + (-self)

    def __mul__(self, other):
        if other == 0:
            return 0
        elif other == 1:
            return self
        elif not isinstance(other, PolynomialNTTRepresentation):
            raise NotImplementedError(TYPE_MISMATCH_ERR)
        elif self.modulus != other.modulus:
            raise NotImplementedError(MODULUS_MISMATCH_ERR)
        elif self.degree != other.degree or len(self.values) != len(other.values) or len(self.values) != self.degree or len(other.values) != other.degree:
            raise NotImplementedError(DEGREE_MISMATCH_ERR)
        elif self.root_order != other.root_order or self.root_order != 2*self.degree or other.root_order != 2*other.degree:
            raise NotImplementedError(ROOT_ORDER_MISMATCH_ERR)
        elif self.root != other.root:
            raise NotImplementedError(ROOT_MISMATCH_ERR)
        return PolynomialNTTRepresentation(
            modulus=self.modulus,
            degree=self.degree,
            root=self.root,
            inv_root=self.inv_root,
            root_order=self.root_order,
            values=[
                _cent(val=x * y, mod=self.modulus, halfmod=self.halfmod, logmod=self.logmod)
                for x, y in zip(self.values, other.values)],)

    def __rmul__(self, other):
        return self.__mul__(other=other)


def transform(
    x: Union[PolynomialCoefficientRepresentation, PolynomialNTTRepresentation],
) -> Union[PolynomialNTTRepresentation, PolynomialCoefficientRepresentation]:
    if isinstance(x, PolynomialCoefficientRepresentation):
        x_coefs: List[int] = deepcopy(x.coefficients)
        root_powers = [pow(x.root, i, x.modulus) for i in range(x.degree)]
        bit_rev_root_powers = _bit_reverse_copy(val=root_powers)
        _cooley_tukey_ntt(val=x_coefs, mod=x.modulus, halfmod=x.halfmod, logmod=x.logmod, deg=x.degree, root=x.root,
                          brv_powers=bit_rev_root_powers)
        return PolynomialNTTRepresentation(
            modulus=x.modulus,
            degree=x.degree,
            root=x.root,
            inv_root=x.inv_root,
            root_order=x.root_order,
            values=x_coefs,
        )
    elif isinstance(x, PolynomialNTTRepresentation):
        x_vals: List[int] = deepcopy(x.values)
        inv_root_powers = [pow(x.inv_root, i, x.modulus) for i in range(x.degree)]
        bit_rev_inv_root_powers = _bit_reverse_copy(val=inv_root_powers)
        _gentleman_sande_intt(val=x_vals, mod=x.modulus, halfmod=x.halfmod, logmod=x.logmod, deg=x.degree, inv_root=x.inv_root,
                              brv_powers=bit_rev_inv_root_powers)
        return PolynomialCoefficientRepresentation(
            modulus=x.modulus,
            degree=x.degree,
            root=x.root,
            inv_root=x.inv_root,
            root_order=x.root_order,
            coefficients=x_vals,
        )
    else:
        raise NotImplementedError(NTT_NOT_IMPLEMENTED_ERR)

