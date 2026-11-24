import mpmath, json
from mpmath import mp, mpf, sqrt, pi, zeta, log, nstr, catalan
mp.dps = 100
ds = ['-4.42642639776412573896646426700982294597089843','6.44737916742048409591156799003400888601121855',
      '-8.64562722577499762278519416620919538044465012','11.2033347447617946274563894200398608793390327',
      '-14.2544419508879307022434578497774873126461318','17.9099574613079499493133409135313285657740977',
      '-22.2756604016684282054623534553349123908541596','27.2106839535892584638834790532553431850258772']
L = log(2)
cand = [mpf(1), pi, pi**2, pi**3, L, L**2, L**3, catalan, zeta(3), pi*L, pi**2*L, pi*catalan, sqrt(2), sqrt(2)*pi, sqrt(2)*L]
names = ['1','pi','pi^2','pi^3','L','L^2','L^3','C','z3','piL','pi^2L','piC','sq2','sq2pi','sq2L']
for i,dsv in enumerate(ds[:4]):
    T = mpf(dsv)
    rel = mpmath.pslq(cand+[T], maxcoeff=10**8)
    if rel and rel[-1]!=0:
        txt = " + ".join(f"({r})*{nm}" for r,nm in zip(rel,names) if r) + f" + ({rel[-1]})*d{i+1}"
        mc = max(abs(x) for x in rel)
        print(f"d{i+1}: maxc {nstr(mc,3)}: {txt[:160]}")
    else:
        print(f"d{i+1}: NONE")
