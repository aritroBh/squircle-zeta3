
integ2 = '''
import sympy as sp
import sym_G, time
t0=time.time()
G21, G12, syms = sym_G.build()
s0,s1,c0,c1,S0,Sc = syms['s0'],syms['s1'],syms['c0'],syms['c1'],syms['S0'],syms['Sc']
sf,cf,Sh,Ch = sp.symbols('sf cf Sh Ch')
subsd = {
  sp.sqrt(1+Sc**2/S0**2): 1/cf,
  s0: sf*Ch+cf*Sh, c0: cf*Ch-sf*Sh,
  s1: Sh*cf-Ch*sf, c1: -(Ch*cf+Sh*sf),
  S0: cf/Sh, Sc: -sf/Sh,
}
gsum = sp.expand(G21+G12).subs(subsd)
gsum = sp.expand(gsum)
print("prep time", time.time()-t0, "terms:", len(gsum.as_ordered_terms()))
# Reduce sf^2 -> 1-cf^2 to canonical (so only sf^0 or sf^1 remain); even part = sf^0 terms
gsum_red = sp.expand(gsum.subs(sf**2, 1-cf**2))
# iterate terms
terms = sp.Add.make_args(gsum_red)
print("reduced terms:", len(terms))
# Each term: coeff * sf^a * cf^b * Sh^p * Ch^q (a in {0,1}); keep a=0
phi = sp.Symbol('phi')
def integ_sincos(a, b):
    # integral_{-h}^{h} sin^a phi cos^b phi dphi ; a even (0), b integer >=? 
    return sp.integrate(sp.cos(phi)**b, (phi, ))  # a=0
# collect
hh = sp.Symbol('h', positive=True)
tot = 0
nodd=0; nkept=0
monolist=[]
for t in terms:
    d = t.as_powers_dict()
    a = d.get(sf,0); b = d.get(cf,0); p = d.get(Sh,0); q = d.get(Ch,0)
    rest = t / (sf**a * cf**b * Sh**p * Ch**q)
    rest = sp.simplify(rest)
    if a==1:
        nodd+=1
        continue
    nkept+=1
    monolist.append((b,p,q,rest))
print("kept", nkept, "odd-dropped", nodd)
import collections
bp = collections.Counter((b) for b,p,q,r in monolist)
print("cf powers b present:", dict(bp))
pp = collections.Counter((p) for b,p,q,r in monolist)
print("Sh powers p present:", dict(pp))
qq = collections.Counter((q) for b,p,q,r in monolist)
print("Ch powers q present:", dict(qq))
# Save monolist
import pickle
with open('monolist.pkl','wb') as f: pickle.dump([(b,p,q,sp.srepr(r)) for b,p,q,r in monolist],f)
print("saved monolist.pkl; total terms kept", len(monolist))
'''
with open('step_07_integ2.py','w') as f: f.write(integ2)
import subprocess,time
r = subprocess.run(['./.venv/bin/python','step_07_integ2.py'],capture_output=True,text=True,timeout=600)
print(r.stdout); print(r.stderr[-2000:] if r.stderr else "")
