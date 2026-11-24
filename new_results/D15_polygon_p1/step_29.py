import mpmath as mp
mp.mp.dps = 170

def make_coeffs(m):
    A = mp.pi/m
    def coeffs(th):
        ls = [mp.log(mp.sin(th + i*A)) for i in range(m)]
        sv = [mp.sin(th + i*A) for i in range(m)]
        cv = [mp.cos(th + i*A) for i in range(m)]
        T = [mp.fsum([sv[i]*ls[i]**k for i in range(m)]) for k in range(4)]
        Tp = [mp.fsum(cv)] + [mp.fsum([cv[i]*(ls[i]**k + k*ls[i]**(k-1)) for i in range(m)]) for k in range(1,4)]
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
        return [D0, D0*(rho1 + E1), D0*(rho2 + rho1*E1 + E2), D0*(rho3 + rho2*E1 + rho1*E2 + E3)]
    return coeffs

def P_coeffs(m):
    A = mp.pi/m; x = A/2
    coeffs = make_coeffs(m)
    return [4*m*mp.quad(lambda th: coeffs(th)[k], [0, x], method='tanh-sinh') for k in range(4)]

out = {}
for m in [3,4]:
    cm = P_coeffs(m)
    out[m] = [mp.factorial(k)*cm[k] for k in range(4)]
    print(f"m={m}: P'' = {mp.nstr(out[m][2],140)}")
    print(f"      P'''= {mp.nstr(out[m][3],140)}")
import json
with open('vals_m34_170.json','w') as f:
    json.dump({str(m):[mp.nstr(v,170) for v in out[m]] for m in [3,4]}, f)
print("saved")
