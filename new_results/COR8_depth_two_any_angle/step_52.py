
import pickle, re
from fractions import Fraction as Fr
import sympy as sp
c = sp.Symbol('c')
D = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/full_sols_raw.pkl','rb'))
sols, alln = D['sols'], D['alln']

def parse_prod(nm):
    """'Li2((c + 1)**2/(4*c))*log(c)*log(1-c)' -> [('Li2','(c + 1)**2/(4*c)'),('log','c'),('log','1-c')]
       'Li2(x)**2' -> [('Li2','x'),('Li2','x')]"""
    parts = []
    i = 0
    while i < len(nm):
        m = re.match(r'(Li[23]|log)\(', nm[i:])
        kind = m.group(1)
        j = i + len(kind) + 1
        depth = 1
        k = j
        while depth:
            if nm[k] == '(': depth += 1
            elif nm[k] == ')': depth -= 1
            k += 1
        arg = nm[j:k-1]
        parts.append((kind, arg))
        k += 1  # skip ')'
        if nm[k-1:k+3] == '**2':
            parts.append((kind, arg))
            k += 3
        if k < len(nm):
            assert nm[k] == '*'
            k += 1
        i = k
    return parts

terms = []
for j, (kind, nm) in enumerate(alln):
    cf = {E: sols[E][j] for E in (-1,0,1) if sols[E][j]}
    if not cf: continue
    if kind == 'Li4':
        terms.append(('Li4', sp.sympify(nm), cf))
    else:
        parts = [(k, sp.sympify(a)) for k, a in parse_prod(nm)]
        terms.append(('prod', parts, cf))
pickle.dump(terms, open('identity_terms.pkl','wb'))
print(len(terms), "terms stored")
for t in terms[:6]: print(t)
