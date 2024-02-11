# Error Messages

# check the following error messages for redundancy and consistency
_MUST_BE_INT_ERR: str = "must be int"
_MUST_BE_LIST_ERR: str = "must be list"
_MUST_BE_LIST_W_POW_2_LEN_ERR: str = "must be list with length a power of 2"
_MUST_BE_POS_INT_ERR: str = "must be positive int"
_MUST_BE_INT_GEQ_3_ERR: str = "must be int >= 3"
_MUST_BE_HALF_FLOORED_ERR: str = "must be floored half"
_MUST_BE_BIT_LEN_ERR: str = "must be bit length"
_MUST_HAVE_PROU_ERR: str = "modulus must have primitive root of unity"
_NO_PROU_FOUND_ERR: str = "no primitive root of unity found"
_MUST_BE_BOOL_ERR: str = "must be bool"
_MUST_BE_ODD_PRIME_ERR: str = "must be odd prime"
_MUST_BE_CORRECT_ROOT_ERR: str = "must be correct root"
_MUST_BE_CORRECT_INVERSE_ROOT_ERR: str = "must be correct inverse root"
_DEGREE_MISMATCH_ERR: str = "degree mismatch"
_MUST_CONSTRUCT_BRV_POWERS_CORRECTLY_ERR: str = "must construct bit reverse powers correctly"
_INVALID_NTT_INPUT_ERR: str = "invalid ntt input"
_TYPE_MISMATCH_ERR: str = "type mismatch"
_MODULUS_MISMATCH_ERR: str = "modulus mismatch"
_ROOT_MISMATCH_ERR: str = "root mismatch"
_ROOT_ORDER_MISMATCH_ERR: str = "root order mismatch"
_NORM_NOT_IMPLEMENTED_ERR: str = "norm not implemented"
_NTT_NOT_IMPLEMENTED_ERR: str = "ntt not implemented"
_WEIGHT_NOT_IMPLEMENTED_ERR: str = "weight not implemented"
_MUL_BASE_NOT_IMPLEMENTED_ERR: str = "multiplication of base class not implemented"
_INV_ROOT_MISMATCH_ERR: str = "inverse root mismatch"
_INVALID_REP_TYPE_ERR: str = "invalid representation type"
_MUST_BE_NONEMPTY_ERR: str = "must be nonempty"
_DIMENSION_MISMATCH_ERR: str = "dimension mismatch"
_COEFS_NORM_WEIGHT_METHOD_MISSING_ERR: str = "matrix elements must have a coefs_norm_weight method"
_MUST_BE_TUPLE_ERR: str = "must be tuple"
_MUST_BE_ALGEBRAIC_CLASS_ERR: str = "must be algebraic class"
_DECODING_ERR: str = "decoding error"

class _ValidationError(Exception):
    """Base class for validation errors."""
    def __init__(self, message: str):
        self.message = message
        super().__init__(self.message)


class _MustBeIntError(_ValidationError):
    """Raised when a value must be an integer."""
    pass


class _MustBeListError(_ValidationError):
    """Raised when a value must be a list."""
    pass


class _MustBeListWithPow2LenError(_ValidationError):
    """Raised when a list must have a length that is a power of 2."""
    pass


class _MustBePositiveIntError(_ValidationError):
    """Raised when a value must be a positive integer."""
    pass


class _MustBePosPowerOfTwoError(_ValidationError):
    """Raised when a value must be a positive power of 2."""
    pass


class _MustBeIntGEQ3Error(_ValidationError):
    """Raised when a value must be an integer greater than or equal to 3."""
    pass


class _MustBeFlooredHalfError(_ValidationError):
    """Raised when a value must be a floored half."""
    pass


class _MustBeBitLengthError(_ValidationError):
    """Raised when a value must represent a bit length."""
    pass


class _MustHavePROUError(_ValidationError):
    """Raised when modulus must have primitive root of unity."""
    pass


class _NoPROUFoundError(_ValidationError):
    """Raised when no primitive root of unity is found."""
    pass


class _MustBeBoolError(_ValidationError):
    """Raised when a value must be a boolean."""
    pass


class _MustBeOddPrimeError(_ValidationError):
    """Raised when a value must be an odd prime."""
    pass


class _MustBeCorrectRootError(_ValidationError):
    """Raised when a value must be the correct root."""
    pass


class _MustBeCorrectInverseRootError(_ValidationError):
    """Raised when a value must be the correct inverse root."""
    pass


class _DegreeMismatchError(_ValidationError):
    """Raised when there is a degree mismatch."""
    pass


class _MustConstructBRVPowersCorrectlyError(_ValidationError):
    """Raised when bit-reverse vector powers are constructed incorrectly."""
    pass


class _InvalidNTTInputError(_ValidationError):
    """Raised when an input to NTT is invalid."""
    pass


class _TypeMismatchError(_ValidationError):
    """Raised when there is a type mismatch."""
    pass


class _ModulusMismatchError(_ValidationError):
    """Raised when there is a modulus mismatch."""
    pass


class _RootMismatchError(_ValidationError):
    """Raised when there is a root mismatch."""
    pass


class _RootOrderMismatchError(_ValidationError):
    """Raised when there is a root order mismatch."""
    pass


class _NormNotImplementedError(_ValidationError):
    """Raised when a norm computation is not implemented."""
    pass


class _NTTNotImplementedError(_ValidationError):
    """Raised when an NTT computation is not implemented."""
    pass


class _WeightNotImplementedError(_ValidationError):
    """Raised when a weight computation is not implemented."""
    pass


class _MulBaseNotImplementedError(_ValidationError):
    """Raised when multiplication of a base class is not implemented."""
    pass


class _InvRootMismatchError(_ValidationError):
    """Raised when there is an inverse root mismatch."""
    pass


class _InvalidRepTypeError(_ValidationError):
    """Raised when an invalid representation type is encountered."""
    pass


class _MustBeNonemptyError(_ValidationError):
    """Raised when a collection must be nonempty."""
    pass


class _DimensionMismatchError(_ValidationError):
    """Raised when there is a dimension mismatch."""
    pass


class _CoefsNormWeightMethodMissingError(_ValidationError):
    """Raised when matrix elements must have a 'coefs_norm_weight' method."""
    pass


class _MustBeTupleError(_ValidationError):
    """Raised when a value must be a tuple."""
    pass


class _MustBeAlgebraicClassError(_ValidationError):
    """Raised when a value must be an algebraic class."""
    pass


class _DecodingError(_ValidationError):
    """Raised when there is a decoding error."""
    pass