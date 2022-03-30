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


class EmbeddedIssue(BertiniRealError):
    """Indicates that there's an issue with a decomposition being embedded in another"""
    pass

class PleaseRegather(BertiniRealError):
    """
    Indicates that something changed in the code, and you need to re-gather from br's output
    """
    pass
