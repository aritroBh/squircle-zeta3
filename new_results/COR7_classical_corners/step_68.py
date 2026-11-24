import sys, time
sys.path.insert(0, '.')
from fractions import Fraction as Fr
from cor7_pt import QuadField, QuadPoint

t0 = time.time()
# test on the previously-analyzed sqrt17 points + the sqrt2/sqrt5 corners
QF17 = QuadField(17, verbose=True)
tests = [
    (17, "sqrt17-4", (Fr(-4), Fr(1))),
    (17, "(sqrt17-3)/2", (Fr(-3,2), Fr(1,2))),
    (17, "(sqrt17-1)/8", (Fr(-1,8), Fr(1,8))),
    (5, "1/phi", (Fr(-1,2), Fr(1,2))),
    (5, "sqrt5-2 = 1/phi^3", (Fr(-2), Fr(1))),
    (2, "sqrt2-1", (Fr(-1), Fr(1))),
    (2, "(sqrt2-1)^2", (Fr(3), Fr(-2))),
]
QFs = {17: QF17}
for D, nm, c0 in tests:
    if D not in QFs:
        QFs[D] = QuadField(D, verbose=True)
    P = QuadPoint(QFs[D], f"Q(sqrt{D}): {nm}", c0, primes=[2])
    r, labs, vecs = P.letter_span_rank()
    print(f"Q(sqrt{D}) {nm}: letter-span rank {r}  ({time.time()-t0:.1f}s)")
    if vecs:
        print("   ", vecs)
