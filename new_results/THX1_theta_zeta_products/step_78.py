
try:
    import gmpy2
    print("gmpy2 version:", gmpy2.version())
except Exception as e:
    print("no gmpy2:", e)
import flint
print([m for m in dir(flint.fmpq) if 'rec' in m.lower() or 'recon' in m.lower()])
print([m for m in dir(flint) if 'recon' in m.lower()])
import mpmath as mp
print("mpmath backend:", mp.__version__)
