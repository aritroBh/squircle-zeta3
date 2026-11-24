from fractions import Fraction as Fr
# correct Bareiss elimination (cor9_rank.rank_int) vs cor7_pt's buggy one, on the 5 corners + controls
import sys
WD = '<project>/research/08_theory/COR9_golden_corner_closed_forms/'
sys.path.insert(0, WD)
from cor9_pt import CFQuadField, QuadPoint
from cor9_rank import rank_int

def buggy_rank(vecs):
    M = [list(map(Fr, v)) for v in vecs]
    R = 0
    for c in range(len(M[0])):
        p_ = next((r for r in range(R, len(M)) if M[r][c]), None)
        if p_ is None: continue
        M[R], M[p_] = M[p_], M[R]
        q = M[R][c]; M[R] = [x/q for x in M[R]]
        for r in range(len(M)):
            if r != R and M[r][c]:
                M[r] = [a - b*M[R][c] for a, b in zip(M[r], M[R])]
        R += 1
    return R

pts = [(2, "sqrt2-1", (Fr(-1), Fr(1))), (2, "(sqrt2-1)^2", (Fr(3), Fr(-2))), (2, "1/sqrt2", (Fr(0), Fr(1,2))),
       (5, "1/phi", (Fr(-1,2), Fr(1,2))), (5, "1/phi^3", (Fr(-2), Fr(1))),
       (17, "sqrt17-4", (Fr(-4), Fr(1))), (5, "1/phi^2", (Fr(3,2), Fr(-1,2))), (2, "(sqrt2-1)^4", (Fr(17), Fr(-12)))]
pr = {2: [2], 5: [2], 17: [2]}
QFs = {}
for D, nm, c0 in pts:
    if D not in QFs: QFs[D] = CFQuadField(D)
    primes = pr[D] if D != 5 or nm != "1/phi^2" else [2, 5]
    if nm == "(sqrt2-1)^4": primes = [2, 3]
    P = QuadPoint(QFs[D], nm, c0, primes=primes)
    rr, labs, vecs = P.letter_span_rank()
    M = [list(v) for _, v in vecs]
    print(f"{nm:16s}: cor7_pt (buggy) = {rr}, rank_int (correct) = {rank_int(M)}")
open(WD+'progress.log','a').write("cor7_pt letter_span_rank elimination bug documented; 5 corners unaffected (rank 2 both ways).\n")
