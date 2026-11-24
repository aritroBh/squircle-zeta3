import sys, time
sys.path.insert(0, '.')
from fractions import Fraction as Fr
from cor7_pt import QuadField, QuadPoint
from cor6_pipe import W4TERMS
from cor6_ft import ft_delta22

t0 = time.time()
QFs = {}
tests = []
# split-field candidates from the (vacuous) LV filter
for (D, a, b, e) in [(17,-4,1,1),(17,-3,1,2),(17,5,-1,2),(17,-3,1,4),(17,-1,1,4),(17,5,-1,4),
                     (17,-1,1,8),(17,1,1,8),(17,7,-1,16),(17,7,1,16),(41,-5,1,2),(41,-5,1,8)]:
    tests.append((D, f"c0=({a}{'+' if b>=0 else ''}{b}s{D})/{e}", (Fr(a,e), Fr(b,e)), [2]))
# extra dim-3 test points: rational-prime extras in sqrt2/sqrt5
tests += [
    (2, "(sqrt2-1)^3", (Fr(-7), Fr(5)), [2,7]),
    (2, "(sqrt2-1)^4 = 17-12sqrt2", (Fr(17), Fr(-12)), [2,3]),
    (5, "1/phi^2 = (3-sqrt5)/2", (Fr(3,2), Fr(-1,2)), [2,5]),
    (5, "1/(2phi) = (sqrt5-1)/4 (n=10)", (Fr(-1,4), Fr(1,4)), [2,5]),
    (3, "sqrt3/2 (n=3)", (Fr(0), Fr(1,2)), [2,3]),
]
for D, nm, c0, primes in tests:
    if D not in QFs:
        QFs[D] = QuadField(D)
    try:
        P = QuadPoint(QFs[D], f"Q(sqrt{D}) {nm}", c0, primes=primes)
        r, labs, vecs = P.letter_span_rank()
        if r is None:
            print(f"Q(sqrt{D}) {nm}: letter factorization FAILED (not in generated group)")
            continue
        if r <= 2:
            print(f"Q(sqrt{D}) {nm}: rank {r} (!!)")
        else:
            SOm = P.omega_symbol(W4TERMS)
            d22 = ft_delta22(SOm, 2, len(P.gensE))
            print(f"Q(sqrt{D}) {nm}: rank {r}, TRUE d22 {'ZERO' if not d22 else 'NONZERO ('+str(len(d22))+' comp)'}  ({time.time()-t0:.0f}s)", flush=True)
    except Exception as ex:
        print(f"Q(sqrt{D}) {nm}: ERROR {ex}")
