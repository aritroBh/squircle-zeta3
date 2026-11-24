open('progress.log','a').write('step_08: exact kernels k=1..14 at c=1/2,3/5,1/3,4/5 saved (293s total)\n')
step09 = r'''"""step_09.py -- polynomial law: finite differences of kernel residues.
For each key (i,j), d = i+j <= 5, each real pole rho: res(k) for k = d..14.
Test: (d)-th finite difference of res(k) vanishes for k >= d+1 (degree <= d-1 beyond k=d).
Report max |Delta^d res| over the range, and the exact polynomial (Newton series) with spare points.
"""
import sys, pickle
sys.path.insert(0, '.')
from fractions import Fraction as Fr

def qpart(t):   # plain Quad coeff ((an,ad),(bn,bd)) -> (rational a, rational b)
    return Fr(*t[0]), Fr(*t[1])

ANGLES = ['c12', 'c35', 'c13', 'c45']
for name in ANGLES:
    data = pickle.load(open(f'kern_{name}.pkl', 'rb'))
    K = data['kernels']
    print(f"=== {name} (c={Fr(*data['c'])}, r={Fr(*data['r'])}) ===")
    maxbad = Fr(0); ntest = 0; nfail = 0
    for d in range(1, 6):
        keys = [(i, j) for i in range(d, (d-1)//2, -1) for j in [d-i]]
        for key in keys:
            # collect residue sequences for k = d..14
            seqs = {}   # pole -> list of (a,b) for k=d..14
            consts = []
            for k in range(d, 15):
                P = K[k][key]
                for rho, cf in P['real'].items():
                    seqs.setdefault(rho, []).append(qpart(cf))
                consts.append(qpart(P['const']))
            for rho, seq in seqs.items():
                for comp in range(2):   # 0 = rational part, 1 = theta part (should be 0)
                    s = [x[comp] for x in seq]
                    # k ranges d..14; polynomial law for k >= d+1: take k=d+1..14 -> s[1:]
                    sp = s[1:]
                    # d-th finite difference
                    def fdiff(seq, n):
                        for _ in range(n):
                            seq = [seq[i+1]-seq[i] for i in range(len(seq)-1)]
                        return seq
                    dd = fdiff(list(sp), d)
                    ntest += 1
                    if any(x != 0 for x in dd):
                        nfail += 1
                        if comp == 0:
                            print(f"  FAIL key{key} pole{rho} comp{comp}: Delta^{d} = {dd[:4]}")
                    if any(x != 0 for x in fdiff(list(sp), d+1)):
                        pass
            # constant part should be zero
            if any(x != (Fr(0), Fr(0)) for x in consts):
                print(f"  key{key} nonzero const parts: {consts[:3]}")
    print(f"  tested {ntest} (pole,key,comp) sequences: {nfail} failures of Delta^d=0 (k>=d+1, d<=5)")
print("STEP09 DONE")
'''
open('step_09.py','w').write(step09)
import subprocess
PY = '<project>/.venv/bin/python'
r = subprocess.run(['nice','-n','10',PY,'step_09.py'], capture_output=True, text=True)
print(r.stdout)
print(r.stderr[-1000:] if r.returncode else "RC=0")
