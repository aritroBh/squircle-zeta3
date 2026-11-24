import mpmath as mp
mp.mp.dps = 60

def make_coeffs(m):
    A = mp.pi/m
    def s(i, th): return mp.sin(th + i*A)
    def sp(i, th): return mp.cos(th + i*A)
    def coeffs(th):
        # T_k and T_k'
        ls = [mp.log(s(i,th)) for i in range(m)]
        sv = [s(i,th) for i in range(m)]
        cv = [sp(i,th) for i in range(m)]
        T = [mp.fsum([sv[i]*ls[i]**k for i in range(m)]) for k in range(4)]
        Tp = [mp.mpf(0)] + [mp.fsum([cv[i]*(ls[i]**k + k*ls[i]**(k-1)) for i in range(m)]) for k in range(1,4)]
        S0 = T[0]; S0p = Tp[0]
        U1 = T[1]/S0
        U2 = T[2]/(2*S0) - T[1]**2/(2*S0**2)
        U3 = T[3]/(6*S0) - T[1]*T[2]/(2*S0**2) + T[1]**3/(3*S0**3)
        U1p = Tp[1]/S0 - T[1]*S0p/S0**2
        U2p = Tp[2]/(2*S0) - T[2]*S0p/(2*S0**2) - T[1]*Tp[1]/S0**2 + T[1]**2*S0p/S0**3
        U3p = (Tp[3]/(6*S0) - T[3]*S0p/(6*S0**2)
               - (Tp[1]*T[2]/(2*S0**2) + T[1]*Tp[2]/(2*S0**2) - T[1]*T[2]*S0p/S0**3)
               + (T[1]**2*Tp[1]/S0**3 - T[1]**3*S0p/S0**4))
        lS0 = mp.log(S0)
        R1 = lS0 - U1
        R2 = -lS0 + U1 - U2
        R3 = lS0 - U1 + U2 - U3
        Q0 = S0p/S0
        Q1 = Q0 - U1p
        Q2 = -Q0 + U1p - U2p
        Q3 = Q0 - U1p + U2p - U3p
        rho1 = R1
        rho2 = R2 + R1**2/2
        rho3 = R3 + R1*R2 + R1**3/6
        W1 = 2*Q0*Q1; W2 = 2*Q0*Q2 + Q1**2; W3 = 2*Q0*Q3 + 2*Q1*Q2
        B = 1 + Q0**2
        E1 = W1/(2*B)
        E2 = W2/(2*B) - W1**2/(8*B**2)
        E3 = W3/(2*B) - W1*W2/(4*B**2) + W1**3/(16*B**3)
        D0 = (1/S0)*mp.sqrt(B)
        D1 = rho1 + E1
        D2 = rho2 + rho1*E1 + E2
        D3 = rho3 + rho2*E1 + rho1*E2 + E3
        return D0, D0*D1, D0*D2, D0*D3
    return coeffs

def P_coeffs(m):
    A = mp.pi/m; x = A/2
    coeffs = make_coeffs(m)
    res = []
    for k in range(4):
        I = mp.quad(lambda th: coeffs(th)[k], [0, x], method='tanh-sinh')
        res.append(4*m*I)
    return res  # c_k = [d^k] P ; P^{(k)}(1) = k! c_k

print("m=2 validation against R2 (PROVED there):")
c = P_coeffs(2)
r = mp.sqrt(2)
known = [4*r, mp.mpf(0), mp.pi**2*r/3, 6*r*(mp.zeta(3)-mp.pi**2/3)]
for k in range(4):
    deriv = mp.factorial(k)*c[k]
    res = abs(deriv - known[k])
    print(f"  P^({k})(1) = {mp.nstr(deriv,45)}  known={mp.nstr(known[k],40)}  res={mp.nstr(res,4)}")

print("\nm=3:")
c3 = P_coeffs(3)
for k in range(4):
    print(f"  P^({k})(1) = {mp.nstr(mp.factorial(k)*c3[k],45)}")

print("\nm=4:")
c4 = P_coeffs(4)
for k in range(4):
    print(f"  P^({k})(1) = {mp.nstr(mp.factorial(k)*c4[k],45)}")
