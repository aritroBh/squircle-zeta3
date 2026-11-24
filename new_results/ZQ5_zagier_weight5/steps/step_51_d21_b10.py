
import sys; sys.path.insert(0,'.')
import json, time
import mpmath as mp
import w5lib as W
import pipe
D=21
t0=time.time()
C, cands = pipe.phase_enum(D, B=10)
print(f"enum B=10 {time.time()-t0:.0f}s, N={len(cands)}")
t0=time.time()
K = pipe.phase_kernel(D, C, cands)
json.dump(dict(K=K), open(pipe.cand_dir(D)+"/kernel_b10.json","w"))
print(f"kernel {time.time()-t0:.0f}s dim {len(K)}")
import pickle, os
# numerics only for new candidates (reuse vals cache for indices < 78)
old = pipe.load_vals(D, 1150) if os.path.exists(pipe.cand_dir(D)+"/vals_1150.pkl") else {}
t0=time.time()
mp.mp.dps = 1180
newvals = {}
li5 = []
for i, c in enumerate(cands):
    key = f"li5:{i}"
    if i < 78 and key in old:
        v = old[key]
    else:
        v = W.li5_smart(c['x'].emb(1180), 1150)
    li5.append(v); newvals[key] = v
print(f"li5 done {time.time()-t0:.0f}s")
for k in ("log:lE21","log:ls21","log:la","log:l2","pi2","z5","L","T"):
    newvals[k] = old[k]
pickle.dump({k: mp.nstr(v, 1175) for k, v in newvals.items()},
            open(pipe.cand_dir(D)+"/vals_b10.pkl","wb"))
# monkey-load
def load_b10(D2, P=1150):
    raw = pickle.load(open(pipe.cand_dir(D)+"/vals_b10.pkl","rb"))
    oldd = mp.mp.dps; mp.mp.dps = 1250
    out = {k: mp.mpf(v) for k, v in raw.items()}
    mp.mp.dps = oldd
    return out
vals = load_b10(D)
for P in (830, 520):
    t0=time.time()
    R = pipe.phase_lll(D, C, cands, K, vals, P=P, use_atoms=False)
    print(f"P={P}: exact_rels={len(R['exact_rels'])} exact_zeros={len(R['exact_zeros'])} ({time.time()-t0:.0f}s)")
    for mx,a,res in R["exact_rels"][:2]:
        print("   a0=", a[0], " max|a|=", mx, " resid=", mp.nstr(res,2))
