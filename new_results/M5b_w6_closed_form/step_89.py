import os
m5b='<project>/research/08_theory/M5b_w6_closed_form'
script = '''"""W6-N6 final kernel LLL at <digits> (no dependency hunt; deps file optional).
Usage: -B w6_lll_final.py <digits>   -> w6_lll_final_<digits>.json
"""
import json, math, pickle, sys, time
from fractions import Fraction as F
from pathlib import Path
sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import flint
import mpmath as mp
from w4_symbol import Q, R, A, B
from m5b_common import log, qnum
from w6_lll import build, col_values, set_prec

def lll(vals, digits, label):
    n = len(vals)
    C = mp.mpf(10)**digits
    rows = [[int(i == j) for j in range(n)]+[int(mp.nint(C*v))] for i, v in enumerate(vals)]
    log(f'LLL start {label}: rows={n}, digits={digits}, gram=exact')
    t0 = time.monotonic()
    reduced = flint.fmpz_mat(rows).lll(gram='exact')
    dt = time.monotonic()-t0
    relations = []
    for i in range(n):
        vec = [int(reduced[i, j]) for j in range(n)]
        mx = max(map(abs, vec))
        residual = mp.fsum(c*v for c, v in zip(vec, vals))
        if abs(residual) < mp.mpf(10)**(-digits-30) and mx < mp.mpf(10)**(mp.mpf(digits)/n):
            relations.append(dict(vec=vec, maxcoef=mx, residual=mp.nstr(residual, 12)))
    log(f'LLL done {label}: seconds={dt:.2f}; small relations={len(relations)}')
    return relations

def main():
    digits = int(sys.argv[1])
    d = build(digits)
    nk = len(d['kv'])
    kept = list(range(2*nk))
    set_prec(digits+100)
    vals = [mp.mpf(v) for v in col_values(d, kept)]
    relations = lll(vals, digits, f'W6 symbol-kernel enlarged basis @{digits}')
    out = dict(digits=digits, kept=kept, relations=relations)
    (HERE/f'w6_lll_final_{digits}.json').write_text(json.dumps(out, indent=2)+'\\n')
    nrel = [r for r in relations if r['vec'][0] != 0]
    log(f'W6 final LLL @{digits}: target-involving relations={len(nrel)}')
    for r in nrel:
        log(f'  c0={r["vec"][0]} maxcoef={r["maxcoef"]} residual={r["residual"]}')

if __name__ == '__main__':
    main()
'''
open(os.path.join(m5b,'w6_lll_final.py'),'w').write(script)
import py_compile; py_compile.compile(os.path.join(m5b,'w6_lll_final.py'), doraise=True)
print('ok')
# launch 1300 in background
import subprocess, time
PY='<project>/.venv/bin/python'
logf=open(m5b+'/w6_lll_1300.out','w')
p=subprocess.Popen(['nice','-n','10',PY,'-B','w6_lll_final.py','1300'],cwd=m5b,stdout=logf,stderr=subprocess.STDOUT)
open(m5b+'/w6_lll_1300.pid','w').write(str(p.pid))
print('launched pid', p.pid)
