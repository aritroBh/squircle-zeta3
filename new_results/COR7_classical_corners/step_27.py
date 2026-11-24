import os, sys
print(os.getcwd())
print(os.listdir('.'))
BASE = '<project>/research/'
# quick check: python-flint and mpmath available
import mpmath as mp
from flint import fmpz_mat, fmpz
print("mpmath", mp.__version__ if hasattr(mp,'__version__') else 'ok', "| flint ok")
