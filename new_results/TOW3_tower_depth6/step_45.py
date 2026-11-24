work = '<project>/research/08_theory/TOW3_tower_depth6'
prove = r'''"""prove_V56.py -- PROOF of the Newton identities at depths 5,6, NON-ADJACENT keys, all (d,m).
Each identity  e(delta,s) + e(delta,s+1) == Lambda^{(delta)}(s)  is an exact rational identity
in (d,m) over Q(i). Method: clear denominators (product of all Pi-factors), expand, check 0.
Usage: python prove_V56.py DELTA S
"""
import sys, time
import sympy as sp
import tow3lib as T3

d, m = T3.d, T3.m

def identity(delta, s):
    lhs = sp.expand(T3.e(delta, s) + T3.e(delta, s+1))
    rhs = T3.Lam(delta, s)
    t0 = time.time()
    rem = sp.together(lhs - rhs)
    num = sp.fraction(rem)[0]
    num = sp.expand(num)
    is0 = num == 0
    print(f'non-adjacent ({delta},{s}): numerator == 0 : {is0}  ({round(time.time()-t0,1)}s)', flush=True)
    return is0

if __name__ == '__main__':
    delta, s = int(sys.argv[1]), int(sys.argv[2])
    ok = identity(delta, s)
    assert ok
'''
with open(f'{work}/prove_V56.py','w') as f: f.write(prove)
import subprocess, os
# run the two VALUE identities first (the key new results), timing them
for (delta, s) in [(5,0)]:
    t0 = __import__('time').time()
    r = subprocess.run(['nice','-n','10', '<project>/.venv/bin/python', 'prove_V56.py', str(delta), str(s)],
                       cwd=work, capture_output=True, text=True, timeout=1700)
    print(r.stdout, r.stderr[-500:] if r.returncode else '', flush=True)
