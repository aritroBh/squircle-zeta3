
validator = '''
import mpmath as mp
mp.mp.dps = 60
import sympy as sp
from poly_coeffs import regular_zeros

# load symbolic G21,G12
import sym_G
G21, G12, syms = sym_G.build()
order = [syms[k] for k in ['s0','s1','c0','c1','tau1','tau2','tau3','gam1','gam2','gam3','S0','Sc','lS0']]
f21 = sp.lambdify(order, G21, 'mpmath')
f12 = sp.lambdify(order, G12, 'mpmath')

def symvals(m, th):
    zeros = regular_zeros(m)
    A = mp.pi/m
    sr = [mp.sin(th-b) for b in zeros]
    sv = [abs(x) for x in sr]
    cv = [mp.cos(th-b)*mp.sign(x) for b,x in zip(zeros,sr)]
    lv = [mp.log(s) for s in sv]
    N = range(1, m-1)   # non-vanishing indices
    s0 = sv[0]; c0v = cv[0]
    s1 = sv[m-1]; c1v = cv[m-1]
    tau = [mp.fsum([sv[i]*lv[i]**k for i in N]) for k in (1,2,3)]
    gam = [mp.fsum([cv[i]*(lv[i]**k + k*lv[i]**(k-1)) for i in N]) for k in (1,2,3)]
    S0 = mp.fsum(sv); Sc = mp.fsum(cv); lS0 = mp.log(S0)
    return [s0,s1,c0v,c1v,tau[0],tau[1],tau[2],gam[0],gam[1],gam[2],S0,Sc,lS0]

# step_09 extractor (perturbative)
from poly_coeffs import make_coeffs
def make_pert_factory(zeros, w, lams, m):
    def make_pert(sx, sy):
        def coeffs(th):
            sraw = [mp.sin(th - b) for b in zeros]
            sv = [abs(x) for x in sraw]
            cv = [mp.cos(th - b)*mp.sign(x) for b, x in zip(zeros, sraw)]
            ls = [mp.log(s) + lam for s, lam in zip(sv, lams)]
            ls[0] += sx; ls[m-1] += sy
            T = [mp.fsum([w[i]*sv[i]*ls[i]**k for i in range(m)]) for k in range(4)]
            Tp = [mp.fsum([w[i]*cv[i] for i in range(m)])] + \\
                 [mp.fsum([w[i]*cv[i]*(ls[i]**k + k*ls[i]**(k-1)) for i in range(m)]) for k in range(1,4)]
            S0 = T[0]; S0p = Tp[0]
            U1 = T[1]/S0
            U2 = T[2]/(2*S0) - T[1]**2/(2*S0**2)
            U3 = T[3]/(6*S0) - T[1]*T[2]/(2*S0**2) + T[1]**3/(3*S0**3)
            U1p = Tp[1]/S0 - T[1]*S0p/S0**2
            U2p = Tp[2]/(2*S0) - T[2]*S0p/(2*S0**2) - T[1]*Tp[1]/S0**2 + T[1]**2*S0p/S0**3
            U3p = (Tp[3]/(6*S0) - T[3]*S0p/(6*S0**2)
                   - (Tp[1]*T[2]/(2*S0**2) + T[1]*Tp[2]/(2*S0**2) - T[1]*T[2]*S0p/S0**3)
                   + (T[1]**2*Tp[1]/S0**3 - T[1]**3*S0p/S0**4))
            lS0 = mp.log(S0); L = S0p/S0
            R1 = lS0 - U1; R2 = -lS0 + U1 - U2; R3 = lS0 - U1 + U2 - U3
            Q0 = -L
            Q1 = L - U1p; Q2 = -L + U1p - U2p; Q3 = L - U1p + U2p - U3p
            rho1 = R1; rho2 = R2 + R1**2/2; rho3 = R3 + R1*R2 + R1**3/6
            W1 = 2*Q0*Q1; W2 = 2*Q0*Q2 + Q1**2; W3 = 2*Q0*Q3 + 2*Q1*Q2
            B = 1 + Q0**2
            E1 = W1/(2*B); E2 = W2/(2*B) - W1**2/(8*B**2); E3 = W3/(2*B) - W1*W2/(4*B**2) + W1**3/(16*B**3)
            D0 = (1/S0)*mp.sqrt(B)
            return D0*(rho3 + rho2*E1 + rho1*E2 + E3)
        return coeffs
    return make_pert
PTS = [(0,0),(1,0),(0,1),(2,0),(1,1),(0,2),(3,0),(2,1),(1,2),(0,3)]
TERMS = [(i,j) for i in range(4) for j in range(4) if i+j<=3]
MX = mp.matrix([[mp.mpf(sx)**i*mp.mpf(sy)**j for (i,j) in TERMS] for sx,sy in PTS])
def Gs_at(make_pert, th):
    vals = mp.matrix([make_pert(mp.mpf(sx), mp.mpf(sy))(th) for sx,sy in PTS])
    cf = mp.lu_solve(MX, vals)
    d = dict(zip(TERMS, list(cf)))
    return d[(2,1)], d[(1,2)]

print("dps now:", mp.mp.dps)
maxr21 = mp.mpf(0); maxr12 = mp.mpf(0)
import itertools
testpts = [(2,'0.3'),(2,'0.9'),(3,'0.5'),(4,'0.2'),(5,'0.6')]
for m, ths in testpts:
    th = mp.mpf(ths)
    mp_ = make_pert_factory(regular_zeros(m),[1]*m,[0]*m,m)
    g21n, g12n = Gs_at(mp_, th)
    sv = symvals(m, th)
    g21s = f21(*sv); g12s = f12(*sv)
    r21 = abs(g21n-g21s); r12 = abs(g12n-g12s)
    maxr21 = max(maxr21, r21); maxr12 = max(maxr12, r12)
    print(f"m={m} th={ths}: |num-sym| G21={mp.nstr(r21,3)}  G12={mp.nstr(r12,3)}")
print("MAX residuals: G21", mp.nstr(maxr21,3), " G12", mp.nstr(maxr12,3), " at dps", mp.mp.dps)
'''
with open('step_01_validate_G.py','w') as f:
    f.write(validator)
import subprocess
r = subprocess.run(['./.venv/bin/python','step_01_validate_G.py'],capture_output=True,text=True)
print(r.stdout); print(r.stderr[-3000:] if r.stderr else "")
