# Silviana Amethyst
# Spring 2022


"""
Defines derived exception types for bertini_real
"""

class BertiniRealError(Exception):
    """Base class for other exceptions in bertini_real"""
    pass


class SurfaceNotSampled(BertiniRealError):
    """Raised when a surface should be sampled before running the executed code"""
    pass
