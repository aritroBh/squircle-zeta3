"""
H1: the n = 2 case of the cross-polytope formula must reproduce the (independently derived,
planner-verified) diamond coefficients a_k of P/sqrt2 = sum a_k d^k from
research/08_theory/D11_rhombus_resurgence/step1_exact.txt.  Here P/sqrt2 = 4 * S_2/S_2(1).
"""
import json, re
import sympy as sp
D11 = 'research/08_theory/D11_rhombus_resurgence/step1_exact.txt'
mine = json.load(open('routeA_n2_K12.json'))['2']['d']
zsub = {sp.Symbol(f'z{2*j}'): sp.zeta(2*j) for j in range(1, 20)}
for line in open(D11):
    m = re.match(r'a_(\d+) = (.*)', line.strip())
    if not m:
        continue
    k = int(m.group(1))
    if k >= len(mine):
        continue
    theirs = sp.sympify(m.group(2).replace('Z', 'z'))
    ours = 4 * sp.sympify(mine[k]).subs(zsub)
    print(f'a_{k}: ours - D11 = {sp.simplify(sp.expand(ours - theirs))}')
