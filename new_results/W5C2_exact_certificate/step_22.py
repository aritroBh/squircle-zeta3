import flint
print("fmpz methods:", [m for m in dir(flint.fmpz) if 'recon' in m.lower() or 'gcd' in m.lower()])
print("nmod_mat methods:", [m for m in dir(flint.nmod_mat) if not m.startswith('_')])
print("fmpq:", [m for m in dir(flint.fmpq) if 'recon' in m.lower()])
# does fmpz_mat.solve exist with algorithm kwarg?
import inspect
try:
    print(inspect.signature(flint.fmpz_mat.solve))
except Exception as e:
    print(e)
# rational reconstruction test
try:
    from flint import fmpz
    print(hasattr(fmpz(3), 'rational_reconstruction'))
except Exception as e: print(e)