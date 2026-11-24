import sys
print(sys.executable)
import flint, mpmath
print('flint', flint.__version__ if hasattr(flint,'__version__') else 'ok', 'mpmath', mpmath.__version__)
base='research/08_theory/A3c_claude_elementary/Q'
print(open(base+'/verify_quotient.py').read())
