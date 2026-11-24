
import time
try:
    import flint
    print("flint version:", flint.__version__)
    print([m for m in dir(flint.nmod_mat) if not m.startswith('_')])
except Exception as e:
    print("flint error:", e)
