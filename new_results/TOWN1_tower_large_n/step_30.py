import os, sys, pickle
from fractions import Fraction as Fr
base = '<project>'
workdir = os.path.join(base, 'research/08_theory/TOWN1_tower_large_n')
os.chdir(workdir)
sys.path.insert(0, workdir)
import driver, trunc_engine as te
from fastred import QSI

H = pickle.load(open('data/towers_D8_K48.pkl','rb'))
H = {k:v for k,v in H.items() if k != '_meta'}
unq4 = driver.unq4

def binom_poly(x, j):
    # x an integer (possibly negative), j>=0: binom(x,j) = x(x-1)...(x-j+1)/j!
    r = Fr(1)
    for i in range(j):
        r *= Fr(x-i)
    for i in range(1, j+1):
        r /= i
    return r

def poly_eval_newton(ks, ys, k):
    k0 = ks[0]
    cur = list(ys); cols=[cur]
    for j in range(1,len(ks)):
        cur=[cur[t+1]-cur[t] for t in range(len(cur)-1)]; cols.append(cur)
    x = k - k0
    out = cols[0][0]
    for j in range(1, len(ks)):
        out = out + cols[j][0].scale(binom_poly(x, j))
    return out

# Diagonal shift
nshift=0; ntot=0
for p in range(1, 24):
    n=2*p; key=(p,p)
    for delta in [0,1,2,3,4]:
        if delta>2*n-4: continue
        kd=H[key].get(delta,{})
        ksd=sorted(k for k in kd if k>=n+1)
        P=delta//2
        if len(ksd)<P+1: continue
        ys=[unq4(kd[k]) for k in ksd[:P+1]]
        pred=poly_eval_newton(ksd[:P+1], ys, n)
        actual=unq4(kd[n])
        ntot+=1
        if pred!=actual: nshift+=1
print(f"Diagonal shift: k=n OFF polynomial in {nshift}/{ntot} diagonal cells (expect = ntot)")

# Off-diagonal k=n ON polynomial
n_on=0; n_tot=0; bad=[]
for n in range(3,41):
    for key in te.keys_of_degree(n):
        p,q=key
        if p==q: continue
        for delta in range(0, min(8,2*n-4)+1):
            kd=H[key].get(delta,{})
            ksd=sorted(k for k in kd if k>=n)
            P=delta//2
            if len(ksd)<P+2: continue
            ys=[unq4(kd[k]) for k in ksd]
            is_poly,deg_exact,npts,nchecks=driver.fd_test(ys,P)
            n_tot+=1
            if is_poly: n_on+=1
            else: bad.append((key,delta))
print(f"Off-diagonal: k=n ON polynomial in {n_on}/{n_tot} cells (expect = n_tot)")
print("  bad:", bad[:6])
with open('progress.log','a') as f:
    f.write(f"step_09 diagonal shift confirmed {nshift}/{ntot} off-poly; off-diag k=n on-poly {n_on}/{n_tot}\n")
