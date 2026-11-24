import sys, os, pickle, sympy as sp, time
wd = "<project>/research/08_theory/DEC2_decagon_exact_assembly"
os.chdir(wd); sys.path.insert(0, wd)
from esum_engine import esum, Hf, H2f, k, x, L1, L2, L3, AA, z2
t0 = time.time()
def Msym(n, a, b):
    l = k + n + 1
    if (a,b)==(0,0): return 1/l
    if (a,b)==(1,0): return -1/l**2
    if (a,b)==(0,1): return -Hf(l)/l
    if (a,b)==(2,0): return 2/l**3
    if (a,b)==(1,1): return Hf(l)/l**2 - z2/l + H2f(l)/l
    if (a,b)==(0,2): return (Hf(l)**2 + H2f(l))/l
def Jexpr(n, a, b, c):
    M = Msym(n, a, b)
    if c == 1: return sp.expand(-esum(M/k, K=1))
    if c == 2: return sp.expand(2*esum(Hf(k-1)*M/k, K=2))
    if c == 3: return sp.expand(-3*esum((Hf(k-1)**2 - H2f(k-1))*M/k, K=3))
tuples = [(n,a,b,c) for n in range(4) for c in (1,2,3) for a in range(4-c) for b in range(4-c-a)]
print("tuples:", len(tuples))
Jtab = {}
for tp in tuples:
    Jtab[tp] = Jexpr(*tp)
print("built Jtab:", len(Jtab), "in", time.time()-t0, "s")
# compare with K1d's Jtab_t.pkl symbolically
K1d = pickle.load(open('<project>/research/08_theory/K1d_engine_fixed/Jtab_t.pkl','rb'))
print("K1d keys:", len(K1d))
mism = 0
for key in K1d:
    if key not in Jtab:
        print("missing", key); mism += 1; continue
    d = sp.expand(Jtab[key] - K1d[key])
    if d != 0:
        mism += 1
        if mism <= 5: print("MISMATCH", key, d)
print("Jtab symbolic mismatches:", mism)
pickle.dump(Jtab, open('Jtab_dec2.pkl','wb'))
with open('progress.log','a') as f:
    f.write(f"step_05: exact single-ratio J-table rebuilt by symbolic summation ({len(Jtab)} keys); vs K1d Jtab_t.pkl: {mism} symbolic mismatches\n")
