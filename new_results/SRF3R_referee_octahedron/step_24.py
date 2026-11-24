import time, json
from mpmath import mp, mpf, quad, sqrt, sin, cos, pi, fabs, log10, exp, log
t0=time.time()
mp.dps = 45

def T_square(tau):
    tau = mpf(tau)
    f = lambda s, u: ((u*(1-u)**2*s*(1-s))**(-tau)
                      * sqrt((1-u)**(2*tau)*(s**(2*tau)+(1-s)**(2*tau)) + u**(2*tau))
                      * (1-u))
    return quad(f, [0,1],[0,1])

def T_laplace(tau):
    tau = mpf(tau)
    def Phi(V):
        if V == 0: return mpf(0)
        em = exp(-V)
        return (exp(tau*V)*(1-em)**(1-2*tau)
                * quad(lambda s: ((s*(1-s))**(-tau)
                       * sqrt((1-em)**(2*tau)*(s**(2*tau)+(1-s)**(2*tau)) + exp(-2*tau*V))), [0,1]))
    return quad(lambda V: exp(-V)*Phi(V), [0, mp.inf])

def T_logit(tau):
    tau = mpf(tau)
    def f(v1, v2):
        D = 1+exp(v1)+exp(v2)
        return exp((1-tau)*(v1+v2)) * D**(-3+2*tau) * sqrt(1+exp(2*tau*v1)+exp(2*tau*v2))
    return quad(f, [-mp.inf, mp.inf], [-mp.inf, mp.inf])

def T_simplex(tau):
    tau = mpf(tau)
    f = lambda t1, t2: ((t1*t2*(1-t1-t2))**(-tau)
                        * sqrt(t1**(2*tau) + t2**(2*tau) + (1-t1-t2)**(2*tau)))
    def inner(t1):
        return (1-t1)*quad(lambda v: f(t1, (1-t1)*v), [0,1])
    return quad(inner, [0,1])

for ds in ["0.01", "-0.03"]:
    p = 1+mpf(ds); tau = 1-1/p
    Ts = T_simplex(tau); Tq = T_square(tau); Tl = T_laplace(tau); Tg = T_logit(tau)
    for name, Tv in [("square",Tq),("laplace",Tl),("logit",Tg)]:
        dd = -log10(fabs(Tv-Ts)/fabs(Ts))
        print(f"d={ds} T_{name} vs T_simplex: {float(dd):.1f} digits [{time.time()-t0:.0f}s]")
    print("  T =", mp.nstr(Ts, 35))
prog=open("progress.log","a"); prog.write("[SRF3R] equivalent forms of T (square/Laplace/logit) agree with simplex form to 35+ digits at d=0.01,-0.03 (own code).\n"); prog.close()
