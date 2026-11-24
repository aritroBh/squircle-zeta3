import pickle, sympy as sp
from fractions import Fraction as Fr
import sys
sys.path.insert(0, ".")
from oct8 import *

sols = pickle.load(open("sols235.pkl","rb"))
vecfun = vec_of_rational

# letters check: all Li3 symbols should only involve {2,3,5}
letters = set()
G = []
for z in sols:
    Sz = S_Li3(z, vecfun)
    for t in Sz: letters.update(t)
    G.append((z, Sz))
print("letters appearing in Li3 symbols:", sorted(letters))

# targets
targets = {
 "Lam1": (Fr(-1,3), Fr(1,3)),
 "Lam2": (Fr(-1,9), Fr(3,1)),
 "Lam3": (Fr(1,4), Fr(2,5)),
 "Lam4": (Fr(1,10), Fr(5,2)),
}
TS = {}
for name,(x,y) in targets.items():
    S = S_Li21(x, y, vecfun)
    for t in S: letters.update(t)
    TS[name] = S
print("letters incl targets:", sorted(letters))

Q = Quotient(sorted(letters))
print("dim V =", len(Q.letters), " dim V^3 =", Q.dim3, " rank Decomp =", Q.rank_D,
      " dim Q =", Q.dim3 - Q.rank_D, " (expected", (Q.dim3-len(Q.letters))//3, ")")

gens = [Q.coord(Sz) for z, Sz in G]
GM = sp.Matrix.hstack(*[g.T for g in gens])
print("rank of Li3 span in V^3:", GM.rank(), " ; rank Li3+Decomp:", sp.Matrix.hstack(GM, sp.Matrix.hstack(*[d.T for d in Q.D])).rank())

results = {}
for name, S in TS.items():
    t = Q.coord(S)
    ok, coeffs = Q.in_span(t, gens)
    results[name] = ok
    print(name, "in Li3-span mod Decomp:", ok)
    if ok:
        # report nonzero coefficients
        nz = [(str(sols[i]), c) for i, c in enumerate(coeffs) if c != 0]
        print("   coeffs:", nz)

# combination 33 Lam1 + 24 Lam2 + 99 Lam3 + 60 Lam4
comb = sadd(smul(TS["Lam1"],33), smul(TS["Lam2"],24), smul(TS["Lam3"],99), smul(TS["Lam4"],60))
t = Q.coord(comb)
ok, coeffs = Q.in_span(t, gens)
results["combo 33L1+24L2+99L3+60L4"] = ok
print("E_arc combo in span:", ok)
if ok:
    nz = [(str(sols[i]), c) for i, c in enumerate(coeffs) if c != 0]
    print("   coeffs:", nz)
pickle.dump({"results":results}, open("Q_results.pkl","wb"))
