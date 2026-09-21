"""
Make the tests import the in-repo `bertini_real` package even when a copy is also installed
in site-packages, by putting the repo `python/` directory first on sys.path.
"""
import os
import sys

_PYTHON_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if _PYTHON_DIR not in sys.path:
    sys.path.insert(0, _PYTHON_DIR)
