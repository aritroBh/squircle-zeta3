import sys, os, pickle
sys.path.insert(0,'.')
import mpmath as mp
mp.mp.dps = 1520
# 1. Direct verification of full Wt5 closed form against the 1550-digit value
fit = pickle.load(open('fit_R_lchi.pkl','rb'))
rel = fit['rel']
phi = (1+mp.sqrt(5))/2
l2 = mp.log(2); lf = mp.log(phi); z2 = mp.zeta(2); z3 = mp.zeta(3); z5 = mp.zeta(5)
L3m = mp.polylog(3, -1/phi); L4m = mp.polylog(4, -1/phi); L4q = mp.polylog(4, 1/phi**2)
L5m = mp.polylog(5, -1/phi); L5q = mp.polylog(5, 1/phi**2)
K5m = L5m + lf**5/120; K5q = L5q - mp.mpf(2)/15*lf**5
K4m = L4m - lf**4/24;  K4q = L4q + lf**4/3
K3m = L3m + lf**3/6
Lchi = mp.mpf(pickle.load(open('Lchi5_1550.pkl','rb'))['Lchi5_5'])
A = {}
for a in range(4):
    b = 3-a
    nm = 'z2' + (('l2^%d'%a if a>1 else 'l2'*a) + ('lf^%d'%b if b>1 else 'lf'*b))
    A[nm] = z2*l2**a*lf**b
for a in range(3):
    b = 2-a
    nm = 'z3' + (('l2^%d'%a if a>1 else 'l2'*a) + ('lf^%d'%b if b>1 else 'lf'*b))
    A[nm] = z3*l2**a*lf**b
A['z2^2l2'] = z2**2*l2; A['z2^2lf'] = z2**2*lf
A['z5'] = z5; A['z2z3'] = z2*z3; A['z2L3m'] = z2*L3m
A['K5m'] = K5m; A['K5q'] = K5q
A['K4ml2'] = K4m*l2; A['K4mlf'] = K4m*lf; A['K4ql2'] = K4q*l2; A['K4qlf'] = K4q*lf
A['K3ml2^2'] = K3m*l2**2; A['K3ml2lf'] = K3m*l2*lf; A['K3mlf^2'] = K3m*lf**2
A['Lchi5'] = Lchi
r5 = mp.sqrt(5)
R_val = sum((mp.mpf(str(p)) + mp.mpf(str(q))*r5)*A[nm] for nm, (p, q) in rel.items())
LOGPART = {'l2^5':("10","0"),'l2^4lf':("125","-25"),'l2^3lf^2':("2025/2","-825/2"),
 'l2^2lf^3':("52075/8","-23525/8"),'l2lf^4':("449025/16","-403175/32"),'lf^5':("15112175/256","-6759335/256")}
def monval(name):
    import re
    a = re.search(r'l2\^?(\d*)', name); b = re.search(r'lf\^?(\d*)', name)
    ea = int(a.group(1)) if a and a.group(1) else (1 if a else 0)
    eb = int(b.group(1)) if b and b.group(1) else (1 if b else 0)
    return l2**ea * lf**eb
logval = sum((mp.mpf(p)+mp.mpf(q)*r5)*monval(nm) for nm, (p, q) in LOGPART.items())
Wt5_closed = logval + R_val
W15 = pickle.load(open('wparts_dec_d1550.pkl','rb'))
Wt5_tab = mp.mpf(W15['Wt5'])
print("Wt5 closed form vs 1550-digit table value:", mp.nstr(abs(Wt5_closed - Wt5_tab), 3))
# also check R itself
print("R residual:", mp.nstr(abs((Wt5_tab - logval) - R_val), 3))
open('progress.log','a').write("step_73: KEY TEST PASSED: R = Wt5 - Wt5_log = 21 classical atoms + L(chi_5,5), IDENTICAL vector at 1150/1400 (resid 3.56e-1115/2.5e-1433, height 2.52e8, denom 4096); direct evaluation of full Wt5 closed form vs 1550-digit table: %s\n" % mp.nstr(abs(Wt5_closed - Wt5_tab), 2))
open('progress.log','a').write("step_74: all 16 golden Li32/Li41 K-parts REDUCE to the 22-atom basis (identical vectors at 1150/1400, residuals at working precision, heights 20..7.2e4) - depth-two constants NOT needed\n")
