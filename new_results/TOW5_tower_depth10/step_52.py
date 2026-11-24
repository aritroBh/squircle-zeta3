import sys, time, pickle, os
work = '<project>/research/08_theory/TOW5_tower_depth10'
sys.path.insert(0, work)
import trunc_engine as te
from fastred import QSI, ZERO
from fractions import Fraction as Fr
def q4(q): return (q.p0, q.p1, q.p2, q.p3)

# custom build: descend only to nmin
def build_tower_min(k, D, rho_k, nmin):
    tower = {key: list(vals) for key, vals in rho_k.items()}
    for n in range(k-1, nmin-1, -1):
        dmax_n = min(D, 2*n-4)
        for key in te.keys_of_degree(n):
            i, j = key
            def U(kk, dp, _t=tower):
                v = _t.get(kk)
                if v is None: return ZERO
                return v[dp] if dp < len(v) else ZERO
            tower[key] = [te.Lambda_ij(i, j, delta, U) for delta in range(dmax_n+1)]
    return tower

# time the build at k=24, D=10
S2 = te.stirling_table(26)
t0=time.time()
rho24 = te.rho_level(24, 10, S2)
t1=time.time()
tw = build_tower_min(24, 10, rho24, 7)
t2=time.time()
print(f"k=24: rho {round(t1-t0,1)}s, build down to 7 {round(t2-t1,1)}s")
# spot-check depth-8 tower vs TOWN1 D=8 data at k=24
H8 = pickle.load(open('<project>/research/08_theory/TOWN1_tower_large_n/data/towers_D8_K52.pkl','rb'))
H8.pop('_meta')
ok=0; tot=0
for key in te.keys_of_degree(7):
    for delta in range(0, 9):
        got = q4(tw[key][delta]); exp = tuple(H8[key][delta][24])
        tot+=1
        if got==exp: ok+=1
print(f"sanity tower k=24 degree 7 depth<=8 vs TOWN1: {ok}/{tot}")
