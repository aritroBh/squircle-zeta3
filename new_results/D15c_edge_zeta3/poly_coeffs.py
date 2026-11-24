import mpmath as mp

def make_coeffs(zeros, weights, lams):
    def coeffs(th):
        sraw = [mp.sin(th - b) for b in zeros]
        sv = [abs(x) for x in sraw]
        cv = [mp.cos(th - b)*mp.sign(x) for b, x in zip(zeros, sraw)]
        ls = [mp.log(s) + lam for s, lam in zip(sv, lams)]
        w = weights
        T = [mp.fsum([w[i]*sv[i]*ls[i]**k for i in range(len(w))]) for k in range(4)]
        Tp = [mp.fsum([w[i]*cv[i] for i in range(len(w))])] + \
             [mp.fsum([w[i]*cv[i]*(ls[i]**k + k*ls[i]**(k-1)) for i in range(len(w))]) for k in range(1,4)]
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

def regular_zeros(m):
    A = mp.pi/m
    return [mp.mpf(0)] + [-(i*A) for i in range(1,m)]

def regular_P_coeffs(m, kmax=3):
    A = mp.pi/m
    coeffs = make_coeffs(regular_zeros(m), [1]*m, [0]*m)
    x = A/2
    return [4*m*mp.quad(lambda th: coeffs(th)[k], [0, x], method='tanh-sinh') for k in range(kmax+1)]

def edge_arc_coeffs(m, kmax=3):
    """full edge arc [0, A] (2m arcs make the full boundary)"""
    A = mp.pi/m
    coeffs = make_coeffs(regular_zeros(m), [1]*m, [0]*m)
    return [2*m*mp.quad(lambda th: coeffs(th)[k], [0, A], method='tanh-sinh') for k in range(kmax+1)]

def P_coeffs_general(zeros, weights, lams, kmax=3):
    coeffs = make_coeffs(zeros, weights, lams)
    cuts = sorted(set([mp.mpf(0)] + [(b % mp.pi) for b in zeros] + [mp.pi]))
    cc = [cuts[0]]
    for c in cuts[1:]:
        if abs(c - cc[-1]) > mp.mpf('1e-30'): cc.append(c)
    cuts = cc
    tot = [mp.mpf(0)]*(kmax+1)
    for (a,b) in zip(cuts[:-1], cuts[1:]):
        for k in range(kmax+1):
            tot[k] += mp.quad(lambda th: coeffs(th)[k], [a, b], method='tanh-sinh')
    return [2*v for v in tot]

def P_direct(zeros, weights, lams, p):
    d = p - 1
    def f(th):
        sv = [abs(mp.sin(th-b)) for b in zeros]
        S = mp.fsum([w*s**p*mp.exp(d*lam) for w,s,lam in zip(weights, sv, lams)])
        Sp = mp.fsum([w*mp.exp(d*lam)*p*s**(p-1)*mp.cos(th-b)*mp.sign(mp.sin(th-b)) for w,s,b,lam in zip(weights,sv,zeros,lams)])
        rho = S**(-1/p); drho = -(1/p)*S**(-1/p-1)*Sp
        return mp.sqrt(rho**2 + drho**2)
    cuts = sorted(set([mp.mpf(0)] + [(b % mp.pi) for b in zeros] + [mp.pi]))
    cc=[cuts[0]]
    for c in cuts[1:]:
        if abs(c-cc[-1])>mp.mpf('1e-30'): cc.append(c)
    return 2*mp.quad(f, cc)

def fd_derivs(zeros, weights, lams, h=mp.mpf('0.004'), K=6):
    xs = [k*h for k in range(-K,K+1)]
    ys = [P_direct(zeros, weights, lams, 1+x) for x in xs]
    M = mp.matrix([[x**j for j in range(2*K+1)] for x in xs])
    cf = mp.lu_solve(M, mp.matrix(ys))
    return [cf[0], cf[1], 2*cf[2], 6*cf[3]]
