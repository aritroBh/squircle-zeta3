"""w6 identification step 10: each of the 12 integrals N_{key,r} (dps-250 run, D=200) against the 13-element level-2 weight-6 basis;
then Q-linear relations among the 12 N together with the 13 basis elements by LLL (flint).  Exploratory; printed numbers from this script."""
import sys, os, json
sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from w6lib import *
D = 200
mp.mp.dps = D + 15
C = consts(); B = basis13(C)
N = {k: mp.mpf(v) for k, v in json.load(open(os.path.join(W, "Nints_250.json")))["N"].items()}
mp.mp.dps = D + 15
for name, x in N.items():
    rel, res, tol = pslq_try(x, B, D, maxcoeff=10 ** 12)
    if rel is None:
        print(f"N[{name}]: no relation with basis13 (D={D})")
    else:
        show(f"N[{name}]", B, rel, res, tol, D)
