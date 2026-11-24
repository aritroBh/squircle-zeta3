import json
from fractions import Fraction as F
from math import lcm, gcd
M4E = '<project>/research/08_theory/M4e_w5_lchi'
d = json.load(open(M4E+'/lchi_final_1200.json'))
names=d['names']; ca=[F(c) for c in d['ca']]; cb=[F(c) for c in d['cb']]
aL=F(d['aL']); bL=F(d['bL'])
cert=json.load(open(M4E+'/cobracket.json'))
# denominators: classical part + L + N5
dens=[a.denominator for a in ca+ca]+[aL.denominator,bL.denominator]
for term in cert['correction']:
    cf=term['coefficient']
    dens += [F(cf[0]).denominator, F(cf[1]).denominator]
D=lcm(*dens)
print("lcm of ALL denominators (incl N5):", D)
import sympy as sp
print("factor:", sp.factorint(D))
# integer relation: D * (W5 - N5) = D * sum... -> W5 coefficient = D
ints=[]
for a,b in zip(ca,cb):
    ints += [abs(a*D), abs(b*D)]
ints += [abs(aL*D), abs(bL*D)]
for term in cert['correction']:
    cf=term['coefficient']; ints += [abs(F(cf[0])*D), abs(F(cf[1])*D)]
print("max |integer coefficient| in primitive relation:", max(ints))
print("W5 coefficient (= D):", D)
# term counts
nz=sum(1 for a,b in zip(ca,cb) if a or b)
li5args=[nm for nm,a,b in zip(names,ca,cb) if nm.startswith('Li5') and (a or b)]
print("nonzero classical atoms:", nz, "; Li5 args:", len(li5args))
print("extras among Li5:", [nm for nm in li5args if nm not in [n for n,_,_ in __import__('nice_search').solve.__defaults__[0]] if False else nm for nm in li5args if 'sqrt2/4' in nm or 'sqrt2/8' in nm or '-2+2' in nm])
print("Li5 args:", li5args)
print("max |numerator| (classical):", max(max(abs(a.numerator) for a in ca), max(abs(b.numerator) for b in cb)))
open(M4E+'/progress.log','a').write("step_57: summary stats computed\n")
