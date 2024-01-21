# polynomials.py

from typing import Dict, List, Union, Tuple
from api.ntt_api import ntt_api
from algebra.ntt import cent, find_primitive_root, has_primitive_root_of_unity

cached_halfmods: Dict[int, int] = {}
cached_logmods: Dict[int, int] = {}
cached_inv_root: Dict[Tuple[int, int], int] = {}


def _validate_modulus_and_coefficients(modulus, coefficients):
    if not isinstance(modulus, int):
        raise TypeError("modulus must be an int")
    elif not has_primitive_root_of_unity(modulus=modulus, root_order=2 * len(coefficients)):
        raise ValueError("Modulus does not have primitive root of unity of appropriate order.")
    elif not all(isinstance(x, int) for x in coefficients):
        raise TypeError("All coefficients must be integers.")


def _validate_same_ring(modulus_a, representation_a, modulus_b, representation_b) -> bool:
    if modulus_a != modulus_b:
        raise NotImplementedError("Cannot do arithmetic with polynomials with different moduli")
    elif len(representation_a) != len(representation_b):
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
    representation: List[int]

    def __init__(self, modulus: int, representation: List[int]):
        """
        Initialize a PolynomialRepresentation instance.

        Args:
            modulus (int): The modulus q for arithmetic operations.
            representation (List[int]): The d-list of values.

        Raises:
            TypeError: If the modulus is not an int.
            ValueError: If the modulus - 1 is not divisible by twice the
                        number of coefficients, or there is no primitive
                        root of unity of correct degree.
            ValueError: If representation is not a list of integers.
        """
        _validate_modulus_and_coefficients(modulus, representation)
        self.modulus = modulus
        self.representation = [cent(val=x, modulus=self.modulus, halfmod=self.halfmod, logmod=self.logmod) for x in representation]

    def __str__(self):
        return f"PolynomialRepresentation(modulus={self.modulus}, representation={self.representation})"

    def __repr__(self):
        return self.__str__()

    def __iter__(self):
        return iter(self.representation)

    def __eq__(self, other):
        return isinstance(other, PolynomialRepresentation) and len(self.representation) == len(other.representation) and all((x - y) % self.modulus == 0 for x, y in zip(self, other))

    def __add__(self, other):
        # Component-wise addition modulo the modulus
        if other == 0:
            return self
        elif not isinstance(other, PolynomialRepresentation):
            raise NotImplementedError("Addition with PolynomialRepresentation and non-PolynomialRepresentation is not defined.")
        _validate_same_ring(modulus_a=self.modulus, representation_a=self.representation, modulus_b=other.modulus,
                            representation_b=other.representation)
        return PolynomialRepresentation(modulus=self.modulus, representation=[cent(val=x+y, modulus=self.modulus, halfmod=self.halfmod, logmod=self.logmod) for x, y in zip(self.representation, other.representation)])

    def __radd__(self, other):
        if other == 0:
            return self
        return self + other

    def __neg__(self):
        return PolynomialRepresentation(modulus=self.modulus, representation=[-x for x in self.representation])

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return other + (-self)

    def __mul__(self, other):
        # Context-sensitive, not implemented for base class
        _validate_same_ring(modulus_a=self.modulus, representation_a=self.representation, modulus_b=other.modulus,
                            representation_b=other.representation)

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
        return 2*len(self.representation)

    @property
    def degree(self) -> int:
        return len(self.representation)

    @property
    def root(self) -> int:
        return find_primitive_root(modulus=self.modulus, root_order=2*len(self.representation))

    @property
    def inv_root(self) -> int:
        if (self.degree, self.modulus) not in cached_inv_root:
            cached_inv_root[(self.degree, self.modulus)] = pow(base=self.root, exp=self.modulus-2, mod=self.modulus)
        return cached_inv_root[(self.degree, self.modulus)]


class PolynomialCoefficientRepresentation(PolynomialRepresentation):
    def __init__(self, modulus: int, representation: List[int]):
        super().__init__(modulus=modulus, representation=representation)

    def __str__(self):
        result: str = f"PolynomialCoefficientRepresentation(modulus={self.modulus}, coefficient representation="
        for i, c in enumerate(self.representation):
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
        return PolynomialCoefficientRepresentation(modulus=self.modulus, representation=super().__add__(other).representation)

    def __mul__(self, other):
        # compute product by foiling
        super().__mul__(other)
        if not isinstance(other, PolynomialCoefficientRepresentation):
            raise NotImplementedError("Addition with PolynomialCoefficientRepresentation and non-PolynomialCoefficientRepresentation is not defined.")
        c: List[int] = [0 for _ in range(2 * self.degree)]
        for i, x in enumerate(self.representation):
            for j, y in enumerate(other.representation):
                c[i + j] += x * y
                if abs(c[i+j]) > 3*self.modulus//2:
                    c[i+j] = cent(val=c[i+j], modulus=self.modulus, halfmod=self.halfmod, logmod=self.logmod)
        return PolynomialCoefficientRepresentation(modulus=self.modulus, representation=[cent(val=x - y, modulus=self.modulus, halfmod=self.halfmod, logmod=self.logmod) for x, y in zip(c[:self.degree], c[self.degree:])])

    def norm(self, p: Union[int, str]) -> int:
        if p != "infty":
            raise NotImplementedError(f"norm for p={p} not implemented")
        return max(abs(x) for x in self.representation)

    def weight(self) -> int:
        return sum(1 if x % self.modulus != 0 else 0 for x in self.representation)


class PolynomialNTTRepresentation(PolynomialRepresentation):
    def __init__(self, modulus: int, representation: List[int]):
        super().__init__(modulus=modulus, representation=representation)

    def __str__(self):
        return f"PolynomialNTTRepresentation(modulus={self.modulus}, values={self.representation})"

    def __eq__(self, other):
        if not isinstance(other, PolynomialNTTRepresentation):
            return False
        return super().__eq__(other)

    def __add__(self, other):
        if not isinstance(other, PolynomialNTTRepresentation):
            raise NotImplementedError(
                f"Addition for {type(self)} and {type(other)} not implemented"
            )
        return PolynomialNTTRepresentation(modulus=self.modulus, representation=super().__add__(other).representation)

    def __mul__(self, other):
        # compute Hadamard product
        super().__mul__(other)
        if not isinstance(other, PolynomialRepresentation):
            raise NotImplementedError("Addition with PolynomialRepresentation and non-PolynomialRepresentation is not defined.")
        return PolynomialNTTRepresentation(
            modulus=self.modulus,
            representation=[cent(val=x * y, modulus=self.modulus, halfmod=self.halfmod, logmod=self.logmod) for x, y in zip(self.representation, other.representation)]
        )


def transform_coefficients_to_ntt(x: PolynomialCoefficientRepresentation) -> PolynomialNTTRepresentation:
    y: List[int] = ntt_api(val=x.representation, modulus=x.modulus, inv_flag=False)
    return PolynomialNTTRepresentation(modulus=x.modulus, representation=y)


def transform_ntt_to_coefficients(x: PolynomialNTTRepresentation) -> PolynomialCoefficientRepresentation:
    y: List[int] = ntt_api(val=x.representation, modulus=x.modulus, inv_flag=True)
    return PolynomialCoefficientRepresentation(modulus=x.modulus, representation=y)


def transform(
    x: Union[PolynomialCoefficientRepresentation, PolynomialNTTRepresentation],
) -> Union[PolynomialNTTRepresentation, PolynomialCoefficientRepresentation]:
    if not (isinstance(x, PolynomialCoefficientRepresentation) or isinstance(x, PolynomialNTTRepresentation)):
        raise NotImplementedError(f"Transform for {type(x)} not implemented")
    elif isinstance(x, PolynomialCoefficientRepresentation):
        return transform_coefficients_to_ntt(x=x)
    return transform_ntt_to_coefficients(x=x)
