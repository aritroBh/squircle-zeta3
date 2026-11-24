import sys, time
sys.path.insert(0, '<project>/research/08_theory/COR6_true_point_symbols')
from fractions import Fraction as Fr
from cor6_nf import NF
from cor6_table import Point
from cor6_ft import ft_delta22

CANDS17 = [
    ("sqrt17: c0 = sqrt17-4", (Fr(-4), Fr(1))),
    ("sqrt17: c0 = (sqrt17-3)/2", (Fr(-3,2), Fr(1,2))),
    ("sqrt17: c0 = (5-sqrt17)/2", (Fr(5,2), Fr(-1,2))),
    ("sqrt17: c0 = (sqrt17-3)/4", (Fr(-3,4), Fr(1,4))),
    ("sqrt17: c0 = (sqrt17-1)/4", (Fr(-1,4), Fr(1,4))),
    ("sqrt17: c0 = (5-sqrt17)/4", (Fr(5,4), Fr(-1,4))),
    ("sqrt17: c0 = (sqrt17-1)/8", (Fr(-1,8), Fr(1,8))),
    ("sqrt17: c0 = (sqrt17+1)/8", (Fr(1,8), Fr(1,8))),
    ("sqrt17: c0 = (7-sqrt17)/16", (Fr(7,16), Fr(-1,16))),
    ("sqrt17: c0 = (7+sqrt17)/16", (Fr(7,16), Fr(1,16))),
]
t0 = time.time()
for label, c0 in CANDS17:
    try:
        P = Point(label, [-17,0,1], c0, basis="auto", verbose=False)
        SOm = P.omega_symbol()
        d22 = ft_delta22(SOm, P.F.d, len(P.basis))
        print(f"{label}: dimV {P.dimV}, letter-span {P.letter_span_dim()}, TRUE d22 = {'ZERO' if not d22 else 'NONZERO '+str(d22)}  ({time.time()-t0:.1f}s)", flush=True)
    except Exception as ex:
        print(f"{label}: FAILED: {ex}  ({time.time()-t0:.1f}s)", flush=True)
