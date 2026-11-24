import sys, pickle, time, random
WD = '<project>/research/08_theory/COR9_golden_corner_closed_forms/'
sys.path.insert(0, WD)
from fractions import Fraction as Fr
from cor9_pt import CFQuadField, QuadPoint
from cor9_rank import true_rank
from cor9_lib import factorint
survivors = pickle.load(open(WD+"census_filter_survivors.pkl","rb"))
random.seed(9)
sample = []
for D in sorted(survivors):
    for a, b, e in random.sample(survivors[D], min(3, len(survivors[D]))):
        sample.append((D, a, b, e))
def norm_primes(D, a, b, e):
    ps = {2}
    for n in (a*a - D*b*b, (e-a)*(e-a) - D*b*b, (e+a)*(e+a) - D*b*b, e):
        ps |= set(factorint(abs(n)).keys())
    return sorted(ps)
agree, disagree, skipped = 0, 0, 0
QFs = {}
for D, a, b, e in sample:
    if D not in QFs: QFs[D] = CFQuadField(D)
    primes = norm_primes(D, a, b, e)
    try:
        P = QuadPoint(QFs[D], "x", (Fr(a, e), Fr(b, e)), primes=primes)
        r1, _, _ = P.letter_span_rank()
    except ValueError:
        skipped += 1; continue
    r2, M, pids, rel = true_rank(D, a, b, e)
    if r1 == r2: agree += 1
    else:
        disagree += 1
        print("DISAGREE:", D, a, b, e, r1, r2)
print(f"cross-check on {len(sample)} random survivors: agree {agree}, disagree {disagree}, skipped (non-principal) {skipped}")
open(WD+'progress.log','a').write(f"{time.strftime('%H:%M:%S')} step_72: rank cross-check cor7_pt vs cor9_rank on sample: agree {agree}/{agree+disagree}.\n")
