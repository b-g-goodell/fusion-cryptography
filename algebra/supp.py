# supp.py
# supporting functions
from math import ceil as _ceil
from typing import Dict as _Dict, Tuple as _Tuple, Union as _Union

CACHED_PRIMITIVE_ROOTS: _Dict[_Tuple[int, int], int] = {}
CACHED_IS_ODD_PRIME: _Dict[int, bool] = {}
CACHED_HAS_PRIMITIVE_ROOT_OF_UNITY: _Dict[_Tuple[int, int], bool] = {}
CACHED_IS_ROOT_OF_UNITY: _Dict[_Tuple[int, int, int], bool] = {}
CACHED_IS_PRIMITIVE_ROOT_OF_UNITY: _Dict[_Tuple[int, int, int], bool] = {}
CACHED_FIND_PRIMITIVE_ROOT: _Dict[_Tuple[int, int], _Union[int, None]] = {}


def is_odd_prime(val: int) -> bool:
    """
    Check if x is an odd prime number.
    :param val: Input number.
    :type val: int
    :return b: Boolean indicating whether x is an odd prime.
    :rtype: bool
    """
    if not isinstance(val, int):
        raise TypeError(f"Value must be an integer, not {type(val)}")
    elif val not in CACHED_IS_ODD_PRIME:
        CACHED_IS_ODD_PRIME[val] = val >= 3 and all(val % i != 0 for i in range(2, _ceil(val**0.5) + 1))
    return CACHED_IS_ODD_PRIME[val]


def has_primitive_root_of_unity(modulus: int, root_order: int) -> bool:
    """
    Check if a modulus has a primitive root of unity of a given order.
    :param modulus: Integer modulus
    :type modulus: int
    :param root_order: Order of primitive root of unity
    :type root_order: int
    :return: True if modulus has a primitive root of unity of order root_order, False otherwise
    """
    if not isinstance(modulus, int) or not isinstance(root_order, int):
        raise TypeError(f"Modulus and root order must be integers, not {(type(modulus), type(root_order))}")
    elif (modulus, root_order) not in CACHED_HAS_PRIMITIVE_ROOT_OF_UNITY:
        CACHED_HAS_PRIMITIVE_ROOT_OF_UNITY[(modulus, root_order)] = (
                modulus >= 3
                and root_order >= 2
                and (modulus - 1) % root_order == 0)
    return CACHED_HAS_PRIMITIVE_ROOT_OF_UNITY[(modulus, root_order)]


def is_root_of_unity(purported_root: int, modulus: int, root_order: int) -> bool:
    """
    Check if val is a root of unity of order root_order modulo modulus.
    :param purported_root: Input integer
    :type purported_root: int
    :param modulus: Input modulus
    :type modulus: int
    :param root_order: Order of the root of unity
    :type root_order: int
    :return b: Boolean indicating whether val is a root of unity of order root_order modulo modulus
    :rtype: bool
    """
    if not isinstance(purported_root, int) or not isinstance(modulus, int) or not isinstance(root_order, int):
        raise TypeError("Input must be integers")
    if (purported_root, modulus, root_order) not in CACHED_IS_ROOT_OF_UNITY:
        CACHED_IS_ROOT_OF_UNITY[(purported_root, modulus, root_order)] = (
                modulus >= 3
                and root_order >= 2
                and pow(purported_root, root_order, modulus) == 1)
    return CACHED_IS_ROOT_OF_UNITY[(purported_root, modulus, root_order)]


def is_primitive_root(purported_root: int, modulus: int, root_order: int) -> bool:
    """
    Check if val is a primitive root of order root_order modulo modulus.
    :param purported_root: Input integer
    :type purported_root: int
    :param modulus: Input modulus
    :type modulus: int
    :param root_order: Order of the root of unity
    :type root_order: int
    :return b: Boolean indicating whether val is a primitive root of order root_order modulo modulus
    :rtype: bool
    """
    if not isinstance(purported_root, int) or not isinstance(modulus, int) or not isinstance(root_order, int):
        raise TypeError("Input must be integers")
    if (purported_root, modulus, root_order) not in CACHED_IS_PRIMITIVE_ROOT_OF_UNITY:
        is_rou: bool = is_root_of_unity(purported_root=purported_root, modulus=modulus, root_order=root_order)
        is_prim: bool = is_rou and all(pow(purported_root, i, modulus) != 1 for i in range(1, root_order))
        CACHED_IS_PRIMITIVE_ROOT_OF_UNITY[(purported_root, modulus, root_order)] = is_prim
    return CACHED_IS_PRIMITIVE_ROOT_OF_UNITY[(purported_root, modulus, root_order)]


def find_primitive_root(modulus: int, root_order: int) -> int:
    """
    Find a primitive root of order root_order modulo modulus. Naive loop that first checks 2, then 3, then 4...
    :param modulus: Input modulus
    :type modulus: int
    :param root_order: Order of the root of unity
    :type root_order: int
    :return r: root of unity r such that r**root_order % q == 1 and r**i % q != 1 for all i in range(1, root_order)
    :rtype: int
    """
    if not isinstance(modulus, int) or not isinstance(root_order, int):
        raise TypeError("modulus and root_order must be integers")
    elif not has_primitive_root_of_unity(modulus=modulus, root_order=root_order):
        raise ValueError("modulus does not have primitive root of unity of correct order")
    elif (modulus, root_order) not in CACHED_FIND_PRIMITIVE_ROOT:
        CACHED_FIND_PRIMITIVE_ROOT[(modulus, root_order)] = None
        r: int = 2
        while r < modulus and not is_primitive_root(purported_root=r, modulus=modulus, root_order=root_order):
            r += 1
        if not is_primitive_root(purported_root=r, modulus=modulus, root_order=root_order):
            raise RuntimeError(
                f"No primitive root found with modulus={modulus}, root_order={root_order}, despite that"
                + f"has_primitive_root_of_unity returned True."
            )
        CACHED_FIND_PRIMITIVE_ROOT[(modulus, root_order)] = r
    return CACHED_FIND_PRIMITIVE_ROOT[(modulus, root_order)]
