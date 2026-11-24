"""step_21 (TOW5R): the LEAD top-coefficient table beta_6 (weight 7) from kernel data
vs the zigzag prediction implied by TOW5's proved tower law (d=6 even):
  beta(-a) = (-1)^{j+1} C(6,j) dl6 on every key (j = second index);
  beta(-b) = (-1)^j C(6,j) dl6 for j < 3;  beta(0) = (-1)^j C(6,j) dl6/2 on diagonal (3,3);
  everything else (incl. beta(a), all sqrt2 parts) = 0;  dl6 = sigma_5 A_6/6! = 61/720.
beta(r) = leading coefficient (k^5 coeff) of the quintic residue at pole r = 5th fwd diff / 5!.
"""
import pickle
from fractions import Fraction as Fr
from math import comb, factorial
wts = '<project>/research/08_theory/WTS_weight_k_structure/theory'
K = {k: pickle.load(open(f'{wts}/kern_{k}.pkl','rb')) for k in range(7, 14)}
keys6 = [(6,0),(5,1),(4,2),(3,3)]
A6 = 61
dl6 = Fr(A6, factorial(6))   # sigma_5 = (-1)^3 = -1?  sigma_m = (-1)^floor((m+1)/2): m=5 -> floor(3)=3 -> -1
# NOTE: paper defines delta_d = [x^d](sech+tanh) = sigma_{d-1} A_d/d! ; for d=6: sigma_5 = (-1)^3 = -1
dl6 = -Fr(A6, factorial(6))
print('delta_6 =', dl6, flush=True)
def lead(ys):
    cur = list(ys)
    for _ in range(5): cur = [cur[t+1]-cur[t] for t in range(len(cur)-1)]
    return cur[0]  # constant 5th forward difference = top Newton coefficient (paper's beta)
poles = {(Fr(0),Fr(0)): '0', (Fr(1),Fr(1)): 'a', (Fr(-1),Fr(-1)): '-a', (Fr(1),Fr(-1)): '-b'}
bad = 0; checked = 0
report = []
for (i0, j0) in keys6:
    for pole, name in poles.items():
        ys2 = [K[k][(i0,j0)]['real'].get(pole,(Fr(0),)*4) for k in range(7,14)]
        lc = [lead([y[c] for y in ys2]) for c in range(4)]
        # prediction
        if name == '-a': pred0 = (-1)**(j0+1) * comb(6, j0) * dl6
        elif name == '-b' and j0 < 3: pred0 = (-1)**j0 * comb(6, j0) * dl6
        elif name == '0' and i0 == j0: pred0 = (-1)**j0 * comb(6, j0) * dl6 / 2
        else: pred0 = Fr(0)
        got = lc[0] + 2*lc[1]  # value at sqrt2=2? no: compare components: prediction is rational (component 0), comp1..3 = 0
        okc = (lc[0] == pred0 and lc[1] == 0 and lc[2] == 0 and lc[3] == 0)
        checked += 1
        if not okc: bad += 1
        report.append(((i0,j0), name, lc[0], pred0, okc))
for r in report: print('  beta_%s%s = %s  pred %s  %s' % (r[0], r[1], r[2], r[3], 'OK' if r[4] else 'MISMATCH'), flush=True)
print(f'STEP_21 VERDICT: beta_6 table: {checked} cells, bad = {bad}', flush=True)
assert bad == 0
