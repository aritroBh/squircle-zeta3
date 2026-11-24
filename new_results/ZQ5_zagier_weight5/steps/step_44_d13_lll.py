
import sys; sys.path.insert(0,'.')
import json, time
import mpmath as mp
import w5lib as W
import pipe
D=13
C = W.letters_info(D)
J, cands = pipe.load_cands(D)
K = [list(map(int, r)) for r in json.load(open(pipe.cand_dir(D)+"/kernel.json"))["K"]]
t0=time.time()
vals = pipe.phase_numerics(D, C, cands, P=1150)
print(f"numerics {time.time()-t0:.0f}s")
mp.mp.dps = 60
print("L(chi13,5) =", mp.nstr(vals["L"], 40))
print("T = L/sqrt13 =", mp.nstr(vals["T"], 40))
for ua in (False, True):
    for P in (520, 830):
        R = pipe.phase_lll(D, C, cands, K, vals, P=P, use_atoms=ua)
        print(f"P={P} atoms={ua}: exact_rels={len(R['exact_rels'])} exact_zeros={len(R['exact_zeros'])}")
        for mx,a,res in R["exact_rels"][:2]:
            print("   a0=", a[0], " max|a|=", mx, " resid=", mp.nstr(res,2))
