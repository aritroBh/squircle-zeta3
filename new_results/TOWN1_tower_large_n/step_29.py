import os, sys, pickle
base = '<project>'
workdir = os.path.join(base, 'research/08_theory/TOWN1_tower_large_n')
os.chdir(workdir)
sys.path.insert(0, workdir)
import driver, trunc_engine as te
from fastred import QSI

H = pickle.load(open('data/towers_D8_K48.pkl','rb'))
H = {k:v for k,v in H.items() if k != '_meta'}
unq4 = driver.unq4

# Diagonal shift confirmation: for diagonal keys (p,p), fit polynomial on k=n+1.. and
# evaluate at k=n; the k=n value should NOT lie on the polynomial (it's the exceptional initial value).
from math import comb
def poly_eval_newton(ks, ys, k):
    # Newton forward from ks[0]
    k0 = ks[0]
    cur = list(ys); cols=[cur]
    for j in range(1,len(ks)):
        cur=[cur[t+1]-cur[t] for t in range(len(cur)-1)]; cols.append(cur)
    x = k - k0
    out = cols[0][0]
    for j in range(1, len(ks)):
        out = out + cols[j][0].scale(comb(x, j) if x>=j else __import__('fractions').Fraction(comb(x,j)) )
    return out

print("Diagonal shift check (k=n off polynomial for diagonal keys):")
nshift = 0; ntot = 0
for p in range(1, 24):
    n = 2*p; key=(p,p)
    for delta in [0,1,2,3,4]:
        if delta > 2*n-4: continue
        kd = H[key].get(delta,{})
        ksd = sorted(k for k in kd if k>=n+1)
        P = delta//2
        if len(ksd) < P+1: continue
        ys = [unq4(kd[k]) for k in ksd[:P+1]]
        # value predicted by polynomial at k=n
        pred = poly_eval_newton(ksd[:P+1], ys, n)
        actual = unq4(kd[n])
        ntot += 1
        if pred != actual:
            nshift += 1
print(f"  diagonal cells where k=n value is OFF the polynomial: {nshift}/{ntot} (expect all, confirming shift)")

# And: off-diagonal keys SHOULD have k=n ON the polynomial
print("Off-diagonal k=n ON polynomial check:")
n_on=0; n_tot=0
for n in range(3, 41):
    for key in te.keys_of_degree(n):
        p,q=key
        if p==q: continue
        for delta in [0,1,2,3,4]:
            if delta>2*n-4: continue
            kd=H[key].get(delta,{})
            ksd=sorted(k for k in kd if k>=n)
            P=delta//2
            if len(ksd)<P+2: continue
            ys=[unq4(kd[k]) for k in ksd]
            is_poly,deg_exact,npts,nchecks=driver.fd_test(ys,P)
            n_tot+=1
            if is_poly: n_on+=1
print(f"  off-diagonal cells with k=n ON polynomial: {n_on}/{n_tot} (expect all)")
