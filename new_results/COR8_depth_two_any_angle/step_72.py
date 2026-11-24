
import pickle, re
from fractions import Fraction as Fr
S = pickle.load(open('<project>/research/08_theory/COR8_depth_two_any_angle/full_sols_v2.pkl','rb'))
sols, alln = S['sols'], S['alln']
def parse_prod(nm):
    parts = []
    i = 0
    while i < len(nm):
        m2 = re.match(r'(Li[23]|log)\(', nm[i:])
        kind = m2.group(1)
        j = i + len(kind) + 1
        depth = 1; k = j
        while depth:
            if nm[k] == '(': depth += 1
            elif nm[k] == ')': depth -= 1
            k += 1
        arg = nm[j:k-1]
        parts.append((kind, arg))
        if nm[k:k+3] == '**2':
            parts.append((kind, arg)); k += 3
        if k < len(nm):
            assert nm[k] == '*'
            k += 1
        i = k
    return parts
for j, (kind, nm) in enumerate(alln):
    cf = {E: sols[E][j] for E in (-1,0,1) if sols[E][j]}
    if not cf: continue
    if kind == 'prod':
        print(f"  {cf}  {nm}  -->  {parse_prod(nm)}")
