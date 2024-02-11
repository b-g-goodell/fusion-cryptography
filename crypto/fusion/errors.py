_MUST_BE_MATRIX_ERR: str = "must be matrix"
_MUST_BE_POLY_ERR: str = "must be polynomial"
_ELEM_CLASS_MISMATCH_ERR: str = "elementary class mismatch"
_NORM_TOO_LARGE_ERR: str = "key norm too large"
_WGHT_TOO_LARGE_ERR: str = "key weight too large"
_KEYS_NOT_VALID_ERR: str = "key not valid"
_PARAMS_MISMATCH_ERR: str = "params mismatch"
_CHALL_NOT_VALID_ERR: str = "signature challenge not valid"
_LENGTH_MISMATCH: str = "length mismatch"
_AGG_COEFS_NOT_VALID_ERR: str = "aggregation coefficients not valid"
_MUST_BE_PARAMS_ERR: str = "must be params"


class _ValidationError(Exception):
    """Base class for validation errors."""
    def __init__(self, message: str):
        self.message = message
        super().__init__(self.message)


class _InvalidSecurityParameterError(Exception):
    pass


class _ParameterMismatchError(Exception):
    pass
