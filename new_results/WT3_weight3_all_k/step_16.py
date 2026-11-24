
import os, json
BASE = '<project>'
WORK = os.path.join(BASE, 'research/08_theory/WT3_weight3_all_k')
os.chdir(WORK)
with open('progress.log','a') as fh:
    fh.write("step_01: degree-2 kernels k=3..15 via fastred+buildintegrand; control vs ibp_reduce k=3,4 ALL keys PASS; wt3_kernels.json saved\n")

data = json.load(open('wt3_kernels.json'))
K = data['kernels']
from fractions import Fraction as Fr
import sympy as sp
S2 = sp.sqrt(2)
def qsi_sp(pair):
    return sp.Rational(pair[0]) + sp.Rational(pair[1])*S2
poles_order = [('0','0'), ('1','1'), ('-1','-1'), ('1','-1'), ('-1','1')]
names = {'(0, 0)':'0', '(1, 1)':'a', '(-1, -1)':'-a', '(1, -1)':'b', '(-1, 1)':'-b'}
for k in range(3, 16):
    print(f"--- k={k}")
    for key in ['20','11']:
        terms = {str(tuple(p)): c for p, c in K[str(k)][key]['real']}
        row = []
        for p in poles_order:
            c = terms.get(str((tuple(p))))
            row.append(sp.Rational(c[0]) + sp.Rational(c[1])*S2 if c else sp.Integer(0))
        print(f"  R_{key}: 0:{row[0]}   a:{row[1]}   -a:{row[2]}   b:{row[3]}   -b:{row[4]}")
