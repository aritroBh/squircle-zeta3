import sys, os
print(sys.executable)
try:
    import flint
    print("flint OK", flint.__version__ if hasattr(flint,'__version__') else '')
except Exception as e:
    print("flint FAIL", e)
try:
    import mpmath
    print("mpmath OK", mpmath.__version__)
except Exception as e:
    print("mpmath FAIL", e)
# check mk_deg9.tsv exists
p = "<project>/research/08_theory/D13f_sweep/mk_deg9.tsv"
print("mk_deg9.tsv exists:", os.path.exists(p))
