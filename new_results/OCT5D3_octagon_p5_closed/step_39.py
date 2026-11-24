import mpmath
print("mpmath version:", mpmath.__version__)
import mpmath.libmp.backend as be
print("backend:", be.BACKEND)
try:
    import gmpy2
    print("gmpy2 available:", gmpy2.version())
except ImportError:
    print("no gmpy2")
try:
    import flint
    print("flint:", flint.__version__ if hasattr(flint,'__version__') else 'ok')
except ImportError:
    print("no flint")
