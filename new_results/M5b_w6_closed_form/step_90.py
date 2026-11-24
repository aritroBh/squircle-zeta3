import os
m5b='<project>/research/08_theory/M5b_w6_closed_form'
script = '''"""Assemble R6 (classical part of W6) from the kernel relation; stats; save w6_closed_form.json."""
import json, math, sys
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

def main():
    d1 = json.loads((HERE/'w6_lll_final_1300.json').read_text())
    d2 = json.loads((HERE/'w6_lll_final_1550.json').read_text())
    r1 = [x for x in d1['relations'] if x['vec'][0] != 0]
    r2 = [x for x in d2['relations'] if x['vec'][0] != 0]
    assert len(r1) == 1 and len(r2) == 1, (len(r1), len(r2))
    v1, v2 = r1[0]['vec'], r2[0]['vec']
    same = v1 == v2 or [-x for x in v1] == v2
    print('relation identical at 1300 and 1550 digits:', same)
    print('residuals:', r1[0]['residual'], r2[0]['residual'])
    assert same
    kept = d2['kept']; vec = v2
    c0 = vec[0]
    sol = json.loads((HERE/'w6_symbol_solve.json').read_text())
    names = sol['names']; K = sol['Ksat']
    xp0 = [F(c) for c in sol['xp']]; xq0 = [F(c) for c in sol['xq']]
    nk = len(K)
    ca = xp0[:]; cb = xq0[:]
    for col, coef in zip(kept, vec[1:]):
        if coef == 0:
            continue
        j, parity = divmod(col, 2)
        for i in range(len(names)):
            if K[j][i]:
                if parity == 0:
                    ca[i] -= F(coef)*K[j][i]/c0
                else:
                    cb[i] -= F(coef)*K[j][i]/c0
    terms = {nm: (a, b) for nm, a, b in zip(names, ca, cb) if a or b}
    li6 = sorted([nm for nm in terms if nm.startswith('Li6')])
    li5 = sorted([nm for nm in terms if 'Li5(' in nm])
    li32 = sorted([nm for nm in terms if 'Li32' in nm])
    lchi = [nm for nm in terms if 'Lchi' in nm]
    mx = max(max(abs(a.numerator), abs(b.numerator)) for a, b in terms.values())
    dn = math.lcm(*[a.denominator for a, b in terms.values() for a in (a, b)])
    print(f'R6 classical terms: {len(terms)}; Li6 atoms used: {len(li6)}; Li5 factors: {len(li5)}; Li32 factors: {len(li32)}; Lchi terms: {len(lchi)}')
    print('Li6 args:', li6)
    print(f'max numerator={mx}, lcm denominator={dn}')
    out = dict(names=list(terms), ca=[str(a) for a, b in terms.values()],
               cb=[str(b) for a, b in terms.values()],
               stats=dict(terms=len(terms), li6=len(li6), li5=len(li5), li32=len(li32),
                          lchi=len(lchi), max_num=mx, lcm=dn),
               li6_args=li6, relation_resids=[r1[0]['residual'], r2[0]['residual']])
    (HERE/'w6_closed_form.json').write_text(json.dumps(out, indent=1)+'\\n')
    print('saved w6_closed_form.json')

main()
'''
open(os.path.join(m5b,'w6_assemble.py'),'w').write(script)
import py_compile; py_compile.compile(os.path.join(m5b,'w6_assemble.py'), doraise=True)
print('w6_assemble.py ok')
