
import mpmath as mp, time
mp.mp.dps = 100
from poly_coeffs import regular_zeros
import sym_G
import sympy as sp
# Use the true perturbation extractor (step_09 style)
def make_pert_factory(zeros, w, lams, m):
    def make_pert(sx, sy):
        def coeffs(th):
            sraw = [mp.sin(th - b) for b in zeros]
            sv = [abs(x) for x in sraw]
            cv = [mp.cos(th - b)*mp.sign(x) for b, x in zip(zeros, sraw)]
            ls = [mp.log(s) + lam for s, lam in zip(sv, lams)]
            ls[0] += sx; ls[m-1] += sy
            T = [mp.fsum([w[i]*sv[i]*ls[i]**k for i in range(m)]) for k in range(4)]
            Tp = [mp.fsum([w[i]*cv[i] for i in range(m)])] +                  [mp.fsum([w[i]*cv[i]*(ls[i]**k + k*ls[i]**(k-1)) for i in range(m)]) for k in range(1,4)]
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
    return d[(2,1)] + d[(1,2)]
def Cm_formula(m):
    h = mp.pi/(2*m); s = mp.sin(h)
    return 48*m*s**4*mp.cos(h)*(3*s**2-1)
print("dps:",mp.mp.dps)
for m in [9,10,12]:
    t0=time.time()
    A = mp.pi/m
    mp_ = make_pert_factory(regular_zeros(m),[1]*m,[0]*m,m)
    I = mp.quad(lambda th: Gs_at(mp_, th), [0, A], method='tanh-sinh')
    coeff = 24*m*I
    pred = Cm_formula(m)
    print(f"m={m} ({time.time()-t0:.0f}s): extractor C_m={mp.nstr(coeff,40)}")
    print(f"     residual vs formula = {mp.nstr(abs(coeff-pred),3)}")
