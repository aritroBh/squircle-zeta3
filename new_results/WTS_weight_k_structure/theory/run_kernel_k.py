"""Compute ALL exact reduced kernels (every key, every degree) for one level k with the fast reducer.
usage: python run_kernel_k.py K   -> writes kern_K.pkl  ({(i,j): plain-tuple PF})
plain PF: {'real': {(r0,r1): (c0,c1,c2,c3)}, 'ipos':[(..4..)], 'ineg':[...], 'const':(..4..)}  Fractions p0+p1*sqrt2+p2*i+p3*sqrt2*i
"""
import sys, os, time, pickle
sys.dont_write_bytecode = True
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
from fastred import *
from buildintegrand import build_integrand_pf

def q4(q):
    return (q.p0, q.p1, q.p2, q.p3)

def plain(P):
    return {'real': {(r.p0, r.p1): q4(c) for r, c in P[0].items()},
            'ipos': [q4(e) for e in P[1]], 'ineg': [q4(e) for e in P[2]], 'const': q4(P[3])}

k = int(sys.argv[1])
fn = os.path.join(HERE, f"binteg_{k}.pkl")
if os.path.exists(fn):
    ipf = pickle.load(open(fn, 'rb'))
else:
    ipf = build_integrand_pf(k)
    pickle.dump(ipf, open(fn, 'wb'))
t0 = time.time()
res = reduce_k_fast(k, integrand_pf=ipf)
t1 = time.time()
out = {key: plain(P) for key, P in res.items()}
pickle.dump(out, open(os.path.join(HERE, f"kern_{k}.pkl"), 'wb'))
print(f"k={k}: reduce {t1-t0:.1f}s keys={sorted(out)}", flush=True)
