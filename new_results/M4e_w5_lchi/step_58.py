import json
from fractions import Fraction as F
M4E = '<project>/research/08_theory/M4e_w5_lchi'
d = json.load(open(M4E+'/lchi_final_1200.json'))
names = d['names']; ca=[F(c) for c in d['ca']]; cb=[F(c) for c in d['cb']]
aL=F(d['aL']); bL=F(d['bL'])
# categorize
groups = {'Li5':[], 'Li4':[], 'w32':[], 'z5':[], 'pure':[], 'other':[]}
for nm,a,b in zip(names,ca,cb):
    if not (a or b): continue
    if nm.startswith('Li5'): groups['Li5'].append((nm,a,b))
    elif 'Li4' in nm: groups['Li4'].append((nm,a,b))
    elif nm=='z5': groups['z5'].append((nm,a,b))
    elif ('Li3' in nm or 'z3' in nm) and ('pi2' in nm or 'Li2' in nm): groups['w32'].append((nm,a,b))
    elif 'Li3' in nm or 'z3' in nm or 'Li2' in nm: groups['w32'].append((nm,a,b))
    else: groups['pure'].append((nm,a,b))
for g in groups:
    print(f"== {g}: {len(groups[g])} terms")
    for nm,a,b in groups[g]:
        print(f"   {nm:28s} {a}  {b}")
print("L(chi8,5):", aL, "+", bL, "sqrt2")
print("total nonzero atoms:", sum(len(groups[g]) for g in groups))
