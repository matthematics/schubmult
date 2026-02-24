def pad_tuple(tup, length):
    """Pad a tuple-like with trailing zeros up to ``length``."""
    return (*tup, *(0,) * (length - len(tup)))
