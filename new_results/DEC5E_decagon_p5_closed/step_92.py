import pickle, re
import mpmath as mp
mp.mp.dps = 1520
raw = pickle.load(open('fit_wt5_raw.pkl','rb'))['rel']
r5 = mp.sqrt(5); phi = (1+r5)/2
l2 = mp.log(2); lf = mp.log(phi); z2 = mp.zeta(2); z3 = mp.zeta(3); z5 = mp.zeta(5)
A = {'l2': l2, 'lf': lf, 'z2': z2, 'z3': z3, 'z5': z5,
     'L3m': mp.polylog(3,-1/phi), 'L4m': mp.polylog(4,-1/phi), 'L4q': mp.polylog(4,1/phi**2),
     'L5m': mp.polylog(5,-1/phi), 'L5q': mp.polylog(5,1/phi**2),
     'Lchi5': mp.mpf(pickle.load(open('Lchi5_1550.pkl','rb'))['Lchi5_5'])}
def aval(nm):
    out = mp.mpf(1)
    for m in re.finditer(r'(Lchi5|L3m|L4m|L4q|L5m|L5q|l2|lf|z2|z3|z5)(\^(\d+))?', nm):
        out *= A[m.group(1)]**(int(m.group(3)) if m.group(3) else 1)
    return out
v = sum((mp.mpf(p)+mp.mpf(q)*r5)*aval(nm) for nm,(p,q) in raw.items())
Wt5_tab = mp.mpf(pickle.load(open('wparts_dec_d1550.pkl','rb'))['Wt5'])
print("raw-expanded Wt5 (28 atoms) vs table:", mp.nstr(abs(v-Wt5_tab),3))
# full P5 assembly numeric value for the record
W = pickle.load(open('wparts_dec_d1550.pkl','rb'))
tot = sum(mp.mpf(W[f'Wt{k}']) for k in range(6))
print("P^(5)(1)/L =", mp.nstr(tot, 60))
Lval = mp.sqrt((25-11*r5)/10)
print("P^(5)(1) =", mp.nstr(tot*Lval, 60))
open('progress.log','a').write("step_77: fully-expanded Wt5 (28 atoms, raw basis) vs 1550-digit table: %s\n" % mp.nstr(abs(v-Wt5_tab),2))
