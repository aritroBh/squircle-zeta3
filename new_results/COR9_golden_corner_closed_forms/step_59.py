import sys, pickle, time, json
WD = '<project>/research/08_theory/COR9_golden_corner_closed_forms/'
sys.path.insert(0, WD)
from fractions import Fraction as Fr
from cor9_pt import CFQuadField, QuadPoint
from cor9_lib import factorint
survivors = pickle.load(open(WD+"census_filter_survivors.pkl","rb"))

def norm_primes(D, a, b, e):
    ps = {2}
    for n in (a*a - D*b*b, (e-a)*(e-a) - D*b*b, (e+a)*(e+a) - D*b*b, e):
        ps |= set(factorint(abs(n)).keys())
    return sorted(ps)

t0 = time.time()
rank2 = []
skipped = {}
checked = 0
QFs = {}
for D in sorted(survivors):
    QF = CFQuadField(D)
    for a, b, e in survivors[D]:
        primes = norm_primes(D, a, b, e)
        try:
            P = QuadPoint(QF, f"D{D} ({a},{b})/{e}", (Fr(a, e), Fr(b, e)), primes=primes)
        except ValueError as ex:
            skipped.setdefault(D, 0)
            skipped[D] += 1
            continue
        checked += 1
        r, labs, vecs = P.letter_span_rank()
        if r is not None and r <= 2:
            rank2.append((D, a, b, e, r, [list(map(str, v)) for _, v in vecs], primes))
    print(f"D={D} done ({time.time()-t0:.0f}s), rank2 so far: {len(rank2)}", flush=True)
print(f"\nchecked: {checked}, skipped (non-principal prime): {sum(skipped.values())} in fields {sorted(skipped)}")
print(f"RANK-2 POINTS: {len(rank2)}")
import mpmath as mp
mp.mp.dps = 50
for D, a, b, e, r, vecs, primes in rank2:
    c1 = (mp.mpf(a) + mp.mpf(b)*mp.sqrt(D))/e
    print(f"  Q(sqrt{D}): c0 = ({a}{'+' if b>=0 else ''}{b}s{D})/{e} = {mp.nstr(c1,20)}  alpha = {mp.nstr(2*mp.acos(c1),20)}")
json.dump(rank2, open(WD+"census_rank2_exact.json", "w"))
open(WD+'progress.log','a').write(f"{time.strftime('%H:%M:%S')} step_58: exact rank check; {len(rank2)} rank-2 points, {sum(skipped.values())} skipped.\n")
