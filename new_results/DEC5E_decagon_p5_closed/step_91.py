import pickle
from fractions import Fraction as F
# expand K-atoms into raw atoms: K5m=L5m+lf^5/120, K5q=L5q-(2/15)lf^5,
# K4m=L4m-lf^4/24, K4q=L4q+lf^4/3, K3m=L3m+lf^3/6; K*x*log stays product of raw with logs.
R = pickle.load(open('fit_R_lchi.pkl','rb'))['rel']
LOG = {'l2^5': (F(10),F(0)), 'l2^4lf': (F(125),F(-25)), 'l2^3lf^2': (F(2025,2),F(-825,2)),
       'l2^2lf^3': (F(52075,8),F(-23525,8)), 'l2lf^4': (F(449025,16),F(-403175,32)),
       'lf^5': (F(15112175,256),F(-6759335,256))}
raw = dict(LOG)  # atom -> (p,q)
def add(nm, p, q):
    op, oq = raw.get(nm, (F(0),F(0)))
    raw[nm] = (op+p, oq+q)
    if raw[nm] == (F(0),F(0)): del raw[nm]
for nm, (ps, qs) in R.items():
    p, q = F(ps), F(qs)
    if nm == 'K5m': add('L5m', p, q); add('lf^5', p/120, q/120)
    elif nm == 'K5q': add('L5q', p, q); add('lf^5', -2*p/15, -2*q/15)
    elif nm == 'K4ml2': add('L4ml2', p, q); add('l2lf^4', -p/24, -q/24)
    elif nm == 'K4mlf': add('L4mlf', p, q); add('lf^5', -p/24, -q/24)
    elif nm == 'K4ql2': add('L4ql2', p, q); add('l2lf^4', p/3, q/3)
    elif nm == 'K4qlf': add('L4qlf', p, q); add('lf^5', p/3, q/3)
    elif nm == 'K3ml2^2': add('L3ml2^2', p, q); add('l2^2lf^3', p/6, q/6)
    elif nm == 'K3ml2lf': add('L3ml2lf', p, q); add('l2lf^4', p/6, q/6)
    elif nm == 'K3mlf^2': add('L3mlf^2', p, q); add('lf^5', p/6, q/6)
    else: add(nm, p, q)
print("fully-expanded Wt5, %d atoms:" % len(raw))
for nm in sorted(raw):
    print("  %-12s (%s, %s)" % (nm, raw[nm][0], raw[nm][1]))
import json
pickle.dump({'rel': {k: (str(v[0]), str(v[1])) for k, v in raw.items()}}, open('fit_wt5_raw.pkl','wb'))
