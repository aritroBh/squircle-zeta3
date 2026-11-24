code = '''
import sys, pickle
sys.path.insert(0,'.')
import sympy as sp
from fastred import *
from fractions import Fraction as Fr
v, u = sp.symbols('v u'); I = sp.I; S2 = sp.sqrt(2)
d = sp.Symbol('d', positive=True, integer=True)

# (B2) Int factor for general d: int u^{-2d} du = -u^{-(2d-1)}/(2d-1)
integ = sp.integrate(u**(-2*d), u)
print("(B2) int u^(-2d) du =", sp.simplify(integ), " == -u^(1-2d)/(2d-1):",
      sp.simplify(integ - (-u**(1-2*d)/(2*d-1)))==0)

# (B1) product-top rule, exact check for d=2..12 with FORMAL Laurent series:
# G = sum_{j=1}^{2d-1} g_j u^{-j} (indeterminate coeffs), DL1 actual; then
# [u^{-(2d-2)}](G*DL1) == tau1_0 * g_{2d-1}   (only top of G contributes; needs Lemma A bound j<=2d-1)
tau10 = (-1-I)/2
a = 1+S2; b = S2-1
DL1s = 1/(v-b) + 1/(v+a) - 1/v
okB1 = True
for dd in range(2,13):
    gs = sp.symbols('g1:'+str(2*dd))   # g_1..g_{2d-1}
    G = sum(gs[j-1]*u**(-j) for j in range(1,2*dd))
    prod = sp.expand(G * DL1s.subs(v, I+u))
    ser = sp.series(prod, u, 0, 1).removeO().expand()
    coeff = ser.coeff(u, -(2*dd-2))
    okB1 &= (sp.simplify(coeff - tau10*gs[2*dd-2]) == 0)
print("(B1) [u^{-(2d-2)}](G*DL1) = tau1_0*[u^{-(2d-1)}]G, d=2..12 formal:", okB1)

# (B3) F-branch: F maps -i top to +i top with factor -1 (rev_02b), and DL1 at -i has slope conj(tau1_0)=tau2_0 (rev_02).
# Assemble: F-branch contribution to h_{i,i-1}:  -i * F(Gm_ii DL1)|top+i
#  = -i * ( - (conj(top Gm_ii at order 4i-3) * DL1'(-i)) )  = -i * ( - (-conj(h)/ (4i-3)) * tau2_0 ) = -i*tau2_0*conj(h)/(4i-3)
# B-branch (DL2, same key (i,i)): -i*(Gm_ii DL2)|top+i = -i*(-h/(4i-3))*tau2_0 = i*tau2_0*h/(4i-3)
# total = i*tau2_0*(h - conj h)/(4i-3),  2d-1 = 4i-3 for d=2i-1.  Matches code.
print("(B3) F-branch formula: i*tau2_0*(h_ii - conj h_ii)/(2d-1)  [derived, matches step_39]")

# (V3a/V3b) re-verify recursion + stabilization vs data (independent re-run of step_39)
record = pickle.load(open('record_W2b.pkl','rb'))
tau1 = QSI(-1,0,-1,0)/QSI(2); tau2 = QSI(-1,0,1,0)/QSI(2)
def h(k,key):
    i,j=key; dd=i+j
    rec=record.get((k,key))
    if rec is None: return None
    H=rec['H'][1]; top=2*dd-2
    return H[top-1] if len(H)>=top else ZERO
bad=[]
for k in range(3,13):
    for i in range(1,k+1):
        for j in range(0,i+1):
            dd=i+j
            if dd<2 or k==dd: continue
            hv = h(k,(i,j))
            if hv is None: continue
            s = ZERO
            hA = h(k,(i+1,j))
            if hA is not None: s = s + hA*QSI(i+1)*tau1
            if j==i-1:
                hd = h(k,(i,i))
                if hd is not None: s = s + (hd - hd.conj_i())*QSI(i)*tau2
            elif i>j:
                hB = h(k,(i,j+1))
                if hB is not None: s = s + hB*QSI(j+1)*tau2
            s = s/QSI(2*dd-1)
            if hv != s: bad.append((k,(i,j)))
print("(V3a) top-order recursion h_ij = [(i+1)tau1 h_{i+1,j} + ...]/(2d-1) vs data, k=3..12:", "ALL MATCH" if not bad else bad[:6])
'''
open('rev_03.py','w').write(code)
import sys, subprocess
work = '<project>/research/08_theory/WT2r_review'
r = subprocess.run([sys.executable, f'{work}/rev_03.py'], capture_output=True, text=True, cwd=work)
print(r.stdout); print(r.stderr[-500:] if r.returncode else "")
