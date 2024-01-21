# polynomials.py

from typing import Dict, List, Union, Tuple
from api.ntt_api import ntt_api
from algebra.ntt import cent, find_primitive_root, has_primitive_root_of_unity, is_odd_prime, is_pow_two_geq_two

cached_halfmods: Dict[int, int] = {}
cached_logmods: Dict[int, int] = {}
cached_inv_root: Dict[Tuple[int, int], int] = {}


def _validate_modulus_and_vals(modulus, vals):
    if not isinstance(modulus, int):
        raise TypeError("modulus must be an int")
    elif not isinstance(vals, list) or not all(isinstance(x, int) for x in vals):
        raise TypeError("coefficients must be list of integers")
    elif not is_odd_prime(val=modulus):
        raise ValueError("modulus must be an odd prime")
    elif not has_primitive_root_of_unity(modulus=modulus, root_order=2 * len(vals)):
        raise ValueError("modulus does not have primitive root of unity of appropriate order.")
    elif not is_pow_two_geq_two(val=len(vals)):
        raise ValueError("coefficient vector must have power of 2 length")


def _validate_same_ring(mod_a, vals_a, mod_b, vals_b) -> bool:
    if mod_a != mod_b:
        raise NotImplementedError("Cannot do arithmetic with polynomials with different moduli")
    elif len(vals_a) != len(vals_b):
        raise NotImplementedError("Cannot do arithmetic (addition, multiplication, equality checks) ith polynomials with different max degrees")


class PolynomialRepresentation:
    """
    A base class to represent a polynomial for arithmetic operations modulo (q, X^d + 1).
    This class holds the modulus and the polynomial representation (a list of integers).

    Attributes:
        modulus (int): The modulus q for arithmetic operations.
        representation (List[int]): The representation (d-list of integers))
    """
    modulus: int
    vals: List[int]

    def __init__(self, modulus: int, vals: List[int]):
        """
        Initialize a PolynomialRepresentation instance.

        Args:
            modulus (int): The modulus q for arithmetic operations.
            vals (List[int]): The d-list of values.

        Raises:
            TypeError: If the modulus is not an int.
            ValueError: If the modulus - 1 is not divisible by twice the
                        number of coefficients, or there is no primitive
                        root of unity of correct degree.
            ValueError: If representation is not a list of integers.
        """
        _validate_modulus_and_vals(modulus, vals)
        self.modulus = modulus
        self.vals = [cent(val=x, modulus=self.modulus, halfmod=self.halfmod, logmod=self.logmod) for x in vals]

    def __str__(self):
        return f"PolynomialRepresentation(modulus={self.modulus}, representation={self.vals})"

    def __repr__(self):
        return self.__str__()

    def __iter__(self):
        return iter(self.vals)

    def __eq__(self, other):
        return isinstance(other, PolynomialRepresentation) and len(self.vals) == len(other.vals) and all((x - y) % self.modulus == 0 for x, y in zip(self, other))

    def __add__(self, other):
        # Component-wise addition modulo the modulus
        if other == 0:
            return self
        elif not isinstance(other, PolynomialRepresentation):
            raise NotImplementedError("Addition with PolynomialRepresentation and non-PolynomialRepresentation is not defined.")
        _validate_same_ring(mod_a=self.modulus, vals_a=self.vals, mod_b=other.modulus, vals_b=other.vals)
        return PolynomialRepresentation(modulus=self.modulus, vals=[
            cent(val=x + y, modulus=self.modulus, halfmod=self.halfmod, logmod=self.logmod) for x, y in
            zip(self.vals, other.vals)])

    def __radd__(self, other):
        if other == 0:
            return self
        return self + other

    def __neg__(self):
        return PolynomialRepresentation(modulus=self.modulus, vals=[-x for x in self.vals])

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return other + (-self)

    def __mul__(self, other):
        # Context-sensitive, not implemented for base class
        _validate_same_ring(mod_a=self.modulus, vals_a=self.vals, mod_b=other.modulus, vals_b=other.vals)

    def __rmul__(self, other):
        return self.__mul__(other=other)

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

    @property
    def root_order(self) -> int:
        return 2*len(self.vals)

    @property
    def degree(self) -> int:
        return len(self.vals)

    @property
    def root(self) -> int:
        return find_primitive_root(modulus=self.modulus, root_order=2*len(self.vals))

    @property
    def inv_root(self) -> int:
        if (self.degree, self.modulus) not in cached_inv_root:
            cached_inv_root[(self.degree, self.modulus)] = pow(base=self.root, exp=self.modulus-2, mod=self.modulus)
        return cached_inv_root[(self.degree, self.modulus)]


class PolynomialCoefficientRepresentation(PolynomialRepresentation):
    norm_p: Union[int, str]

    def __init__(self, modulus: int, vals: List[int]):
        super().__init__(modulus=modulus, vals=vals)

    def __str__(self):
        result: str = f"PolynomialCoefficientRepresentation(modulus={self.modulus}, coefficient representation="
        for i, c in enumerate(self.vals):
            if c % self.modulus != 0:
                if i != 0:
                    result += f"{c}X**{i} + "
                else:
                    result += f"{c} + "
        result = result[:-3]
        result +=  ")"
        return result

    def __eq__(self, other):
        if not isinstance(other, PolynomialCoefficientRepresentation):
            return False
        return super().__eq__(other)

    def __add__(self, other):
        if not isinstance(other, PolynomialCoefficientRepresentation):
            raise NotImplementedError(
                f"Addition for {type(self)} and {type(other)} not implemented"
            )
        return PolynomialCoefficientRepresentation(modulus=self.modulus, vals=super().__add__(other).vals)

    def __mul__(self, other):
        # compute product by foiling
        super().__mul__(other)
        if not isinstance(other, PolynomialCoefficientRepresentation):
            raise NotImplementedError("Addition with PolynomialCoefficientRepresentation and non-PolynomialCoefficientRepresentation is not defined.")
        c: List[int] = [0 for _ in range(2 * self.degree)]
        for i, x in enumerate(self.vals):
            for j, y in enumerate(other.vals):
                c[i + j] += x * y
                if abs(c[i+j]) > 3*self.modulus//2:
                    c[i+j] = cent(val=c[i+j], modulus=self.modulus, halfmod=self.halfmod, logmod=self.logmod)
        return PolynomialCoefficientRepresentation(modulus=self.modulus, vals=[
            cent(val=x - y, modulus=self.modulus, halfmod=self.halfmod, logmod=self.logmod) for x, y in
            zip(c[:self.degree], c[self.degree:])])

    @property
    def norm(self) -> int:
        if self.norm_p != "infty":
            raise NotImplementedError(f"norm for p={self.norm_p} not implemented")
        return max(abs(x) for x in self.vals)

    @property
    def weight(self) -> int:
        return sum(1 if x % self.modulus != 0 else 0 for x in self.vals)


class PolynomialNTTRepresentation(PolynomialRepresentation):
    def __init__(self, modulus: int, vals: List[int]):
        super().__init__(modulus=modulus, vals=vals)

    def __str__(self):
        return f"PolynomialNTTRepresentation(modulus={self.modulus}, values={self.vals})"

    def __eq__(self, other):
        if not isinstance(other, PolynomialNTTRepresentation):
            return False
        return super().__eq__(other)

    def __add__(self, other):
        if not isinstance(other, PolynomialNTTRepresentation):
            raise NotImplementedError(
                f"Addition for {type(self)} and {type(other)} not implemented"
            )
        return PolynomialNTTRepresentation(modulus=self.modulus, vals=super().__add__(other).vals)

    def __mul__(self, other):
        # compute Hadamard product
        super().__mul__(other)
        if not isinstance(other, PolynomialRepresentation):
            raise NotImplementedError("Addition with PolynomialRepresentation and non-PolynomialRepresentation is not defined.")
        return PolynomialNTTRepresentation(modulus=self.modulus, vals=[
            cent(val=x * y, modulus=self.modulus, halfmod=self.halfmod, logmod=self.logmod) for x, y in
            zip(self.vals, other.vals)])


def transform_coefficients_to_ntt(x: PolynomialCoefficientRepresentation) -> PolynomialNTTRepresentation:
    transformed_vals: List[int] = ntt_api(val=x.vals, modulus=x.modulus, inv_flag=False)
    return PolynomialNTTRepresentation(modulus=x.modulus, vals=transformed_vals)


def transform_ntt_to_coefficients(x: PolynomialNTTRepresentation) -> PolynomialCoefficientRepresentation:
    transformed_vals: List[int] = ntt_api(val=x.vals, modulus=x.modulus, inv_flag=True)
    return PolynomialCoefficientRepresentation(modulus=x.modulus, vals=transformed_vals)


def transform(x: Union[PolynomialCoefficientRepresentation, PolynomialNTTRepresentation]) -> Union[PolynomialNTTRepresentation, PolynomialCoefficientRepresentation]:
    if not (isinstance(x, PolynomialCoefficientRepresentation) or isinstance(x, PolynomialNTTRepresentation)):
        raise NotImplementedError(f"Transform for {type(x)} not implemented")
    elif isinstance(x, PolynomialCoefficientRepresentation):
        return transform_coefficients_to_ntt(x=x)
    return transform_ntt_to_coefficients(x=x)
