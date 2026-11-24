import json
from fractions import Fraction as F
from math import lcm
import sympy as sp
M4E = '<project>/research/08_theory/M4e_w5_lchi'
d = json.load(open(M4E+'/lchi_final_1200.json'))
names=d['names']; ca=[F(c) for c in d['ca']]; cb=[F(c) for c in d['cb']]
aL=F(d['aL']); bL=F(d['bL'])
cert=json.load(open(M4E+'/cobracket.json'))
dens=[a.denominator for a in ca+ca]+[aL.denominator,bL.denominator]
for term in cert['correction']:
    cf=term['coefficient']; dens += [F(cf[0]).denominator, F(cf[1]).denominator]
D=lcm(*dens)
ints=[]
for a,b in zip(ca,cb): ints += [abs(a*D), abs(b*D)]
ints += [abs(aL*D), abs(bL*D)]
for term in cert['correction']:
    cf=term['coefficient']; ints += [abs(F(cf[0])*D), abs(F(cf[1])*D)]
print("lcm all denominators D:", D, "=", sp.factorint(D))
print("max |integer coefficient|:", max(ints), " W5 coefficient (=D):", D)
nz=sum(1 for a,b in zip(ca,cb) if a or b)
li5args=[nm for nm,a,b in zip(names,ca,cb) if nm.startswith('Li5') and (a or b)]
extra_names={'Li5((-2+2*sqrt2))','Li5((17-12*sqrt2))','Li5((1/2-1/4*sqrt2))','Li5((1/2+1/4*sqrt2))','Li5((1/2-1/2*sqrt2))','Li5((1/2-3/8*sqrt2))'}
s0_used=[nm for nm in li5args if nm not in extra_names]
ex_used=[nm for nm in li5args if nm in extra_names]
print("nonzero classical atoms:", nz, "; Li5 args total:", len(li5args))
print("S0 Li5 used:", s0_used)
print("extras used:", ex_used)
print("max |numerator| (classical):", max(max(abs(a.numerator) for a in ca), max(abs(b.numerator) for b in cb)))
open(M4E+'/progress.log','a').write("step_57: summary stats computed\n")
