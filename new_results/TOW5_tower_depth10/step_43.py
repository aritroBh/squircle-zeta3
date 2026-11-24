import sys, time, pickle, os
work = '<project>/research/08_theory/TOW5_tower_depth10'
sys.path.insert(0, work)
import trunc_engine as te
from fastred import QSI
def q4(q): return (q.p0, q.p1, q.p2, q.p3)
H8 = pickle.load(open('<project>/research/08_theory/TOWN1_tower_large_n/data/towers_D8_K52.pkl','rb'))
H8.pop('_meta')
S2 = te.stirling_table(45)
t0=time.time()
ok = 0; tot = 0
for d in [9, 10, 12]:
    rho_d = te.rho_level(d, 8, S2)
    for key in te.keys_of_degree(d):
        for delta in range(0, 9):
            got = q4(rho_d[key][delta])
            exp = tuple(H8[key][delta][d])
            tot += 1
            if got == exp: ok += 1
print(f"sanity rho depth<=8 vs TOWN1: {ok}/{tot} match ({round(time.time()-t0,1)}s)")
for d in [24, 32, 40]:
    t1=time.time()
    rho_d = te.rho_level(d, 10, S2)
    print(f"rho_level(d={d}, D=10): {round(time.time()-t1,1)}s")
