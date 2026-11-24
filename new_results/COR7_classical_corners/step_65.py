import sys, time
sys.path.insert(0, '.')
from fractions import Fraction as Fr
from cor7_pt import QuadPoint, omega_symbol_pt
from cor6_pipe import W4TERMS
from cor6_ft import ft_delta22

t0 = time.time()
for label, D, c0 in [
    ("sqrt17: c0 = sqrt17-4", 17, (Fr(-4), Fr(1))),
    ("sqrt17: c0 = (sqrt17-3)/2", 17, (Fr(-3,2), Fr(1,2))),
    ("sqrt17: c0 = (sqrt17-1)/8", 17, (Fr(-1,8), Fr(1,8))),
]:
    P = QuadPoint(label, D, c0)
    r, vecs = P.letter_span_rank()
    print(f"{label}: gens {P.gens}, letter-span rank {r}", flush=True)
    print("   letter vecs:", vecs)
    SOm = omega_symbol_pt(P, W4TERMS)
    d22 = ft_delta22(SOm, 2, len(P.gens))
    print("   TRUE d22:", "ZERO" if not d22 else d22, f"({time.time()-t0:.1f}s)")
