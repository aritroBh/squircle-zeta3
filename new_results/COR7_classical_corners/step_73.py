import sys, time
sys.path.insert(0, '.')
from fractions import Fraction as Fr
from cor7_pt import QuadField, QuadPoint
from cor6_pipe import W4TERMS
from cor6_ft import ft_delta22

QF = QuadField(17)
P = QuadPoint(QF, "(7+sqrt17)/16", (Fr(7,16), Fr(1,16)), primes=[2])
r, labs, vecs = P.letter_span_rank()
print("rank:", r, vecs)
SOm = P.omega_symbol(W4TERMS)
d22 = ft_delta22(SOm, 2, len(P.gensE))
print("TRUE d22:", "ZERO" if not d22 else f"NONZERO ({len(d22)} comp)")
