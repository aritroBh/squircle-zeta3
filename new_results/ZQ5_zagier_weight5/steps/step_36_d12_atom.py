
import sys; sys.path.insert(0,'.')
import json, time, pickle
import mpmath as mp
import w5lib as W
import pipe
D=12
C = W.letters_info(D)
J, cands = pipe.load_cands(D)
K = [list(map(int, r)) for r in json.load(open(pipe.cand_dir(D)+"/kernel.json"))["K"]]
vals = pipe.load_vals(D, 1150)
t0=time.time()
R = pipe.phase_lll(D, C, cands, K, vals, P=830, use_atoms=True)
print(f"atom LLL 830: {time.time()-t0:.0f}s rels={len(R['rels'])} exact_rels={len(R['exact_rels'])} exact_zeros={len(R['exact_zeros'])}")
for mx,a,res in R["exact_rels"][:2]:
    print("  exact rel a0=", a[0], " max|a|=", mx, " resid=", mp.nstr(res,2))
t0=time.time()
R5 = pipe.phase_lll(D, C, cands, K, vals, P=520, use_atoms=True)
print(f"atom LLL 520: {time.time()-t0:.0f}s exact_rels={len(R5['exact_rels'])}")
if R5["exact_rels"] and R["exact_rels"]:
    print("two-precision identical:", R5["exact_rels"][0][1] == R["exact_rels"][0][1])
pickle.dump(dict(exact_rels=[(int(m),a,mp.nstr(r,5)) for m,a,r in R["exact_rels"]],
                 exact_zeros=[a for _,a,_ in R["exact_zeros"]], names=R["names"]),
            open(pipe.cand_dir(D)+"/atom_lll.pkl","wb"))
print("saved atom_lll.pkl")
