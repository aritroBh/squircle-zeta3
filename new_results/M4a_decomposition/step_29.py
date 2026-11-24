import subprocess, sys, os, json, datetime
BASE = '<project>'
# verify diff with paper output
r = subprocess.run(['diff', 'ibp_reduce_k4.out', os.path.join(BASE,'paper/squircle_note/outputs/ibp_reduce_k4.out')],
                   capture_output=True, text=True)
print('diff returncode:', r.returncode, '| diff output:', repr(r.stdout[:300]))

code = r'''
import sys, json
sys.dont_write_bytecode = True
sys.path.insert(0, '.')
import sympy as sp
from ibp_reduce import reduce_k, v, S2
from qsqrt2 import QS2
from fractions import Fraction as Fr

a_, b_ = 1 + S2, S2 - 1
letters = {'0': sp.Integer(0), 'a': a_, '-a': -a_, 'b': b_, '-b': -b_}

out = reduce_k(4, verbose=False)

def qs2_of_sympy(e):
    e = sp.expand(e)
    c0 = sp.simplify(e.subs(S2, 0))
    c1 = sp.simplify((e - c0) / S2)
    assert sp.simplify(e - (c0 + c1*S2)) == 0
    return [str(Fr(c0.p, c0.q)), str(Fr(c1.p, c1.q))]

data = {}
for key, R in out.items():
    i, j = key
    if i == 0 and j == 0:
        continue
    terms = []
    for name, z in letters.items():
        c = sp.simplify(R * (v - z))
        c = sp.radsimp(sp.simplify(c.subs(v, z)))
        if c != 0:
            terms.append((name, qs2_of_sympy(c)))
    # verify reconstruction
    recon = sum(sp.Rational(Fr(c[0][0]).numerator, Fr(c[0][0]).denominator)
                + sp.Rational(Fr(c[1][0]).numerator, Fr(c[1][0]).denominator)*S2
                for _, c in [])
    s = sum((sp.Rational(Fr(cq[0]).numerator, Fr(cq[0]).denominator)
             + sp.Rational(Fr(cq[1]).numerator, Fr(cq[1]).denominator)*S2) / (v - letters[name])
            for name, cq in terms)
    assert sp.simplify(s - R) == 0, (key, R, terms)
    data[f'{i},{j}'] = terms

json.dump(data, open('omega4.json','w'), indent=1)
for k in sorted(data, key=lambda s:(-sum(map(int,s.split(","))),s)):
    print(k, data[k])
'''
r = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True, timeout=1200)
print(r.stdout)
print('ERR', r.stderr[:1000])
