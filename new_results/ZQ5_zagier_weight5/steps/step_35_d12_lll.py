
import sys; sys.path.insert(0,'.')
import json, time
import mpmath as mp
import w5lib as W
import pipe
D=12
C = W.letters_info(D)
J, cands = pipe.load_cands(D)
K = [list(map(int, r)) for r in json.load(open(pipe.cand_dir(D)+"/kernel.json"))["K"]]
t0=time.time()
vals = pipe.phase_numerics(D, C, cands, P=1150)
print(f"numerics {time.time()-t0:.0f}s")
for P in (520, 830):
    t0=time.time()
    R = pipe.phase_lll(D, C, cands, K, vals, P=P, use_atoms=False)
    print(f"P={P}: rels={len(R['rels'])} exact_rels={len(R['exact_rels'])} exact_zeros={len(R['exact_zeros'])} ({time.time()-t0:.0f}s)")
    if R["exact_rels"]:
        mx,a,res = R["exact_rels"][0]
        print("   best exact: a0=", a[0], " max|a|=", mx, " resid=", mp.nstr(res,2))
