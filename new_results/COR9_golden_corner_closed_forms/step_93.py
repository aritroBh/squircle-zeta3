import sys, pickle, time, json
WD = '<project>/research/08_theory/COR9_golden_corner_closed_forms/'
sys.path.insert(0, WD)
from fractions import Fraction as Fr
from cor9_pt import CFQuadField, QuadPoint
from cor9_lib import factorint
survivors = pickle.load(open(WD+"census_filter_survivors_v2.pkl","rb"))
def norm_primes(D, a, b, e):
    ps = {2}
    for n in (a*a - D*b*b, (e-a)*(e-a) - D*b*b, (e+a)*(e+a) - D*b*b, e):
        ps |= set(factorint(abs(n)).keys())
    return sorted(ps)
t0 = time.time()
rank2_pt = []
skipped = {}
checked = 0
for D in sorted(survivors):
    QF = CFQuadField(D)
    for a, b, e in survivors[D]:
        primes = norm_primes(D, a, b, e)
        try:
            P = QuadPoint(QF, f"D{D} ({a},{b})/{e}", (Fr(a, e), Fr(b, e)), primes=primes)
        except ValueError:
            skipped[D] = skipped.get(D, 0) + 1
            continue
        checked += 1
        rr, labs, vecs = P.letter_span_rank()
        # NOTE: cor7_pt's elimination is buggy; use exact rank of the (certified) vecs instead
        from cor9_rank import rank_int
        rr2 = rank_int([list(v) for _, v in vecs]) if vecs else rr
        if rr2 is not None and rr2 <= 2:
            rank2_pt.append((D, a, b, e, rr2))
    print(f"D={D} ({time.time()-t0:.0f}s)", flush=True)
print(f"cor7_pt engine: checked {checked}, skipped (non-principal) {sum(skipped.values())}, rank-2 found: {len(rank2_pt)}")
for r in rank2_pt: print("  ", r)
open(WD+'progress.log','a').write(f"{time.strftime('%H:%M:%S')} step_91: cor7_pt engine on 862 (exact rank of certified vecs); rank-2 found {len(rank2_pt)}.\n")
