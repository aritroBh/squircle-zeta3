import os
m5b='<project>/research/08_theory/M5b_w6_closed_form'
script = '''"""Assemble R5'' (classical part of W5'') from the LLL relation; stats; save w5pp_closed_form.json.
Coefficients: atom_i gets particular + kernel lift from the relation; L( chi8,5) separate.
"""
import json, math, sys
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from m5b_common import log

def main():
    d1 = json.loads((HERE/'w5pp_lll_1100.json').read_text())
    d2 = json.loads((HERE/'w5pp_lll_1300.json').read_text())
    r1 = [x for x in d1['relations'] if x['vec'][0] != 0]
    r2 = [x for x in d2['relations'] if x['vec'][0] != 0]
    assert len(r1) == 1 and len(r2) == 1
    v1, v2 = r1[0]['vec'], r2[0]['vec']
    assert d1['kept'] == d2['kept'], 'kept columns differ between precisions'
    same = v1 == v2 or [ -x for x in v1] == v2
    print('relation identical at 1100 and 1300 digits:', same)
    print('residuals:', r1[0]['residual'], r2[0]['residual'])
    assert same
    kept = d2['kept']; vec = v2
    c0 = vec[0]
    sol = json.loads((HERE/'w5pp_symbol_solve.json').read_text())
    names = sol['names']; K = sol['Ksat']
    xp0 = [F(c) for c in sol['xp']]; xq0 = [F(c) for c in sol['xq']]
    nk = len(K)
    ca = xp0[:]; cb = xq0[:]
    aL = F(0); bL = F(0)
    for col, coef in zip(kept, vec[1:]):
        if coef == 0:
            continue
        j, parity = divmod(col, 2)
        if j < nk:
            for i in range(len(names)):
                if K[j][i]:
                    if parity == 0:
                        ca[i] -= F(coef)*K[j][i]/c0
                    else:
                        cb[i] -= F(coef)*K[j][i]/c0
        else:
            if parity == 0:
                aL -= F(coef)/c0
            else:
                bL -= F(coef)/c0
    terms = {nm: (a, b) for nm, a, b in zip(names, ca, cb) if a or b}
    n_li5 = [nm for nm in terms if nm.startswith('Li5')]
    mx = max(max(abs(a.numerator), abs(b.numerator)) for a, b in terms.values())
    dn = math.lcm(*[a.denominator for a, b in terms.values() for a in (a, b)],
                  aL.denominator, bL.denominator)
    print(f\"R5'' classical terms: {len(terms)}; Li5 atoms used: {len(n_li5)}\")
    print('Li5 args:', n_li5)
    print(f'max numerator={mx}, lcm denominator={dn}')
    print(f'L(chi8,5) coefficient: {aL} + {bL} sqrt2')
    out = dict(names=list(terms), ca=[str(a) for a, b in terms.values()],
               cb=[str(b) for a, b in terms.values()], aL=str(aL), bL=str(bL),
               stats=dict(terms=len(terms), li5_used=len(n_li5), max_num=mx, lcm=dn),
               relation_resids=[r1[0]['residual'], r2[0]['residual']])
    (HERE/'w5pp_closed_form.json').write_text(json.dumps(out, indent=1)+'\\n')
    print('saved w5pp_closed_form.json')

main()
'''
open(os.path.join(m5b,'w5pp_assemble.py'),'w').write(script)
import subprocess
PY='<project>/.venv/bin/python'
r=subprocess.run(['nice','-n','10',PY,'-B','w5pp_assemble.py'],capture_output=True,text=True,cwd=m5b,timeout=600)
print(r.stdout); print(r.stderr[-500:])
