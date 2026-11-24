import os, subprocess
d17 = "research/08_theory/D17_zeta5_polygons"
venv_py = "./.venv/bin/python"

script = r'''
import sympy as sp, sys, time, pickle
sys.path.insert(0, "D17DIR")
import sym_D5
t0=time.time()
D, syms = sym_D5.build_D(5, reduced=True)
L0,L1 = syms['L0'], syms['L1']
s0,s1,c0,c1,S0,Sc = syms['s0'],syms['s1'],syms['c0'],syms['c1'],syms['S0'],syms['Sc']
D5 = sp.expand(D[5])
Gs = [sp.expand(D5.coeff(L0,a).coeff(L1,b)) for (a,b) in [(4,1),(3,2),(2,3),(1,4)]]
Gsum = sp.expand(Gs[0]+Gs[1]+Gs[2]+Gs[3])
print("Gsum terms (pre-sub):", len(Gsum.as_ordered_terms()), " build %.0fs"%(time.time()-t0))

sf,cf,Sh,Ch = sp.symbols('sf cf Sh Ch')
subsd = {
  sp.sqrt(1+Sc**2/S0**2): 1/cf,
  s0: sf*Ch+cf*Sh, c0: cf*Ch-sf*Sh,
  s1: Sh*cf-Ch*sf, c1: -(Ch*cf+Sh*sf),
  S0: cf/Sh, Sc: -sf/Sh,
}
t0=time.time()
g = Gsum.subs(subsd)
print("sub done %.0fs"%(time.time()-t0))
g = sp.expand(g)
print("expanded terms:", len(g.as_ordered_terms()), "%.0fs"%(time.time()-t0))
# reduce sf^2 -> 1-cf^2 repeatedly
t0=time.time()
g2 = g
for pw in [2,3,4,5,6,7,8]:
    g2 = sp.expand(g2.subs(sf**pw, sf**(pw-2)*(1-cf**2)))
print("sf-reduced %.0fs"%(time.time()-t0), "terms:", len(g2.as_ordered_terms()))
# classify terms
terms = sp.Add.make_args(g2)
import collections
even=[]; odd=0
for t in terms:
    dd = t.as_powers_dict()
    a = dd.get(sf,0)
    b = dd.get(cf,0); p = dd.get(Sh,0); q = dd.get(Ch,0)
    rest = sp.simplify(t/(sf**a*cf**b*Sh**p*Ch**q))
    if a==0:
        even.append((b,p,q,rest))
    elif a==1:
        odd+=1
    else:
        print("UNEXPECTED sf power", a)
print("even terms:", len(even), " odd(sf^1) terms dropped:", odd)
# check odd part coefficient sums to zero
# (if Gsum even in sf, odd terms should be absent)
bp = collections.Counter(b for b,p,q,r in even)
pp = collections.Counter(p for b,p,q,r in even)
qq = collections.Counter(q for b,p,q,r in even)
print("cf powers b:", dict(bp))
print("Sh powers p:", dict(pp))
print("Ch powers q:", dict(qq))
with open('monolist5.pkl','wb') as fo:
    pickle.dump([(b,p,q,sp.srepr(r)) for b,p,q,r in even], fo)
print("saved monolist5.pkl")
'''.replace("D17DIR", d17)
with open(os.path.join(d17,"subst5.py"),"w") as f: f.write(script)
r = subprocess.run([venv_py, os.path.join(d17,"subst5.py")], capture_output=True, text=True, cwd=d17, timeout=1200)
print(r.stdout); print("ERR:", r.stderr[-3000:])
