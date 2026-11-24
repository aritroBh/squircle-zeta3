import pickle
CWD = "<project>/research/08_theory/G2R_allN_proof_referee"
G2 = "<project>/research/08_theory/G2_M4alpha/"
a = pickle.load(open(CWD + "/run_base/mu4_zsum.pkl", "rb"))
b = pickle.load(open(G2 + "mu4_zsum.pkl", "rb"))
import sympy as sp
same = True
for p in b:
    wa = {tuple(w): sp.simplify(sp.sympify(c)) for w, c in a[p]}
    wb = {tuple(w): sp.simplify(sp.sympify(c)) for w, c in b[p]}
    if set(wa) != set(wb):
        same = False; print("word set mismatch", p)
    for w in wa:
        if sp.simplify(wa[w] - wb[w]) != 0:
            same = False; print("coef mismatch", p, w)
print("regenerated mu4_zsum.pkl == stored artifact:", same)
with open(CWD + "/progress.log", "a") as fh:
    fh.write(f"step_zsum_artifact: run_base/mu4_zsum.pkl == G2 mu4_zsum.pkl: {same}\n")
