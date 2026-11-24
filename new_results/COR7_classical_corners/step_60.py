import sys, time
sys.path.insert(0, '<project>/research/08_theory/COR6_true_point_symbols')
from fractions import Fraction as Fr
from cor6_nf import NF
from cor6_table import Point
from cor6_ft import ft_delta22

CANDS = [
    ("sqrt2: c0 = sqrt2-1", [-2,0,1], (Fr(-1), Fr(1))),
    ("sqrt2: c0 = 3-2sqrt2 = (sqrt2-1)^2", [-2,0,1], (Fr(3), Fr(-2))),
    ("sqrt2: c0 = sqrt2/2 (square, n=4)", [-2,0,1], (Fr(0), Fr(1,2))),
    ("sqrt5: c0 = sqrt5-2 = 1/phi^3", [-5,0,1], (Fr(-2), Fr(1))),
    ("sqrt5: c0 = (sqrt5-1)/2 = 1/phi (golden)", [-5,0,1], (Fr(-1,2), Fr(1,2))),
]
t0 = time.time()
for label, mpoly, c0 in CANDS:
    P = Point(label, mpoly, c0, basis="auto", verbose=False)
    SOm = P.omega_symbol()
    d22 = ft_delta22(SOm, P.F.d, len(P.basis))
    lv = P.letter_vectors()
    print(f"{label}: deg {P.F.d}, dimV {P.dimV}, S0 {P.S0}, letter-span {P.letter_span_dim()}, "
          f"TRUE d22 = {'ZERO' if not d22 else d22}  ({time.time()-t0:.1f}s)", flush=True)
    print("   letter vectors:", {k: tuple(str(x) for x in v) for k, v in lv.items()})
