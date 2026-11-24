import sys, pickle, time, json
WD = '<project>/research/08_theory/COR9_golden_corner_closed_forms/'
sys.path.insert(0, WD)
from fractions import Fraction as Fr
from cor9_pt import CFQuadField, QuadPoint
from cor9_lib import omega_class, ft_delta22, W4TERMS
import cor9_rank; importlib_available = True

rank2 = json.load(open(WD+"census_rank2_complete.json"))
QFs = {2: CFQuadField(2), 5: CFQuadField(5)}
print("=== delta_{2,2} = 0 check for the 5 rank-2 points (exact symbol algebra) ===")
for D, a, b, e, r in rank2:
    P = QuadPoint(QFs[D], f"({a},{b})/{e}", (Fr(a, e), Fr(b, e)), primes=[2])
    rr, labs, vecs = P.letter_span_rank()
    d22 = ft_delta22(P.omega_symbol(W4TERMS), 2, len(P.gensE))
    print(f"  Q(sqrt{D}) c0 = ({a},{b})/{e}: letter rank {rr}, TRUE delta_22 components = {len(d22)}  ->  {'0 ✓' if len(d22)==0 else d22}")

# (rank V, delta) decomposition for the 5 points and the rank-3 survivors
import cor9_rank
import importlib; importlib.reload(cor9_rank)
from cor9_rank import true_rank, letter_vecs, rank_int, kernel_basis
survivors = pickle.load(open(WD+"census_filter_survivors.pkl","rb"))
decomp = {}
for D in sorted(survivors):
    for a, b, e in survivors[D]:
        r, M, pids, rel = true_rank(D, a, b, e)
        rv = rank_int(M)
        delta = r - rv
        # odd-prime part of V
        oddcols = [j for j, pid in enumerate(pids) if pid[0] != 2]
        Modd = [[row[j] for j in oddcols] for row in M]
        rvodd = rank_int(Modd)
        decomp[(r, rv, rvodd, delta)] = decomp.get((r, rv, rvodd, delta), 0) + 1
print("\n(total rank, rank V, rank V_odd, delta) -> count:")
for k in sorted(decomp): print("  ", k, "->", decomp[k])
open(WD+'progress.log','a').write(f"{time.strftime('%H:%M:%S')} step_65: delta_22 checks for 5 points; (rank,V,Vodd,delta) decomposition.\n")
