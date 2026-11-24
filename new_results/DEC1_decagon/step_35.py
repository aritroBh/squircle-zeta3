import os, mpmath as mp, pickle, sympy as sp
os.chdir('research/08_theory/DEC1_decagon')
mp.mp.dps = 100
s5 = mp.sqrt(5)
E3n = pickle.load(open('E3n_num100.pkl','rb'))
def recon(c):
    # c = p + q sqrt(5); pslq [c, 1, s5] -> [r, a, b]: r*c + a + b*s5 = 0
    rel = mp.pslq([c, mp.mpf(1), s5], tol=mp.mpf('1e-90'), maxcoeff=10**9, maxsteps=1000)
    if rel is None: return None
    r,a,b = rel
    return (sp.Rational(-a,r), sp.Rational(-b,r))
exact = {}
maxden = 0
for key,poly in E3n.items():
    row=[]
    for c in poly:
        if abs(c) < mp.mpf('1e-90'):
            row.append((sp.Rational(0),sp.Rational(0))); continue
        rc = recon(c)
        if rc is None:
            print("FAILED recon for", key, mp.nstr(c,20)); row.append(None); continue
        maxden = max(maxden, rc[0].q, rc[1].q)
        row.append(rc)
    exact[key]=row
print("max denominator:", maxden)
# verify reconstruction numerically
worst = mp.mpf(0)
for key,poly in E3n.items():
    for i,c in enumerate(poly):
        p,q = exact[key][i]
        worst = max(worst, abs(c-(mp.mpf(p)+mp.mpf(q)*s5)))
print("reconstruction max residual:", mp.nstr(worst,6))
pickle.dump(exact, open('E3n_exact_pq.pkl','wb'))
# count nonzero
nz = sum(1 for key in exact for i in range(4) if exact[key][i]!=(0,0) and (exact[key][i][0]!=0 or exact[key][i][1]!=0))
print("nonzero t-coefficients:", nz)
with open('progress.log','a') as f:
    f.write(f"step: exact Q(sqrt5) reconstruction of all coeffs; max denom {maxden}, max residual {mp.nstr(worst,4)}\n")