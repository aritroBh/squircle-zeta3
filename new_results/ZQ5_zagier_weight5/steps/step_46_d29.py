
import sys; sys.path.insert(0,'.')
import json, time
import mpmath as mp
import w5lib as W
import pipe
D=29
C, cands = pipe.phase_enum(D, B=8)
pipe.print_cands(D, cands, C)
K = pipe.phase_kernel(D, C, cands)
json.dump(dict(K=K), open(pipe.cand_dir(D)+"/kernel.json","w"))
t0=time.time()
vals = pipe.phase_numerics(D, C, cands, P=1150)
print(f"numerics {time.time()-t0:.0f}s")
for P in (520, 830):
    R = pipe.phase_lll(D, C, cands, K, vals, P=P, use_atoms=False)
    print(f"P={P}: exact_rels={len(R['exact_rels'])} exact_zeros={len(R['exact_zeros'])}")
    for mx,a,res in R["exact_rels"][:2]:
        print("   a0=", a[0], " max|a|=", mx, " resid=", mp.nstr(res,2))
