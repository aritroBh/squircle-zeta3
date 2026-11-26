"""r02 (referee, own code): M_k(alpha) from the paper's definition
   M(u;alpha) = int_0^1 [ (1+g)^(-1-u) sqrt(1 + g^(2-2u) - 2 cos(alpha) g^(1-u)) - 1 ] g^(u-1) dg = sum_k M_k u^k,
by computing, at every quadrature node g, the exact truncated Taylor series in u of the integrand (power-series
arithmetic, log1p/expm1 so that nothing cancels near g=0), and integrating its u^k coefficient with tanh-sinh.
No Cauchy integral, no mu_{k,n}.  usage: python r02_definition.py DPS alpha1 alpha2 ...  (alphas as mpmath expressions)"""
import sys
import mpmath as mp

K = 5

def smul(A, B):
    return [mp.fsum(A[i]*B[k-i] for i in range(k+1)) for k in range(K+1)]

def slog1p_rest(y):         # log(1+y), y[0] == 0
    out = [mp.mpf(0)]*(K+1); P = [mp.mpf(1)] + [mp.mpf(0)]*K
    for t in range(1, K+1):
        P = smul(P, y)
        for k in range(K+1): out[k] += (-1)**(t+1)*P[k]/t
    return out

def sexpm1_rest(y):         # exp(y)-1, y[0] == 0
    out = [mp.mpf(0)]*(K+1); P = [mp.mpf(1)] + [mp.mpf(0)]*K
    for t in range(1, K+1):
        P = smul(P, y)
        for k in range(K+1): out[k] += P[k]/mp.factorial(t)
    return out

def integrand_series(g, cosa):
    Lg = mp.log(g)
    # x(u) = g^(2-2u) - 2 cos(alpha) g^(1-u)
    x = [g**2*(-2*Lg)**k/mp.factorial(k) - 2*cosa*g*(-Lg)**k/mp.factorial(k) for k in range(K+1)]
    x0 = x[0]
    y = [mp.mpf(0)] + [x[k]/(1 + x0) for k in range(1, K+1)]
    half_log = [mp.log1p(x0)/2] + [v/2 for v in slog1p_rest(y)[1:]]
    l1g = mp.log1p(g)
    Lam = list(half_log)
    Lam[0] -= l1g; Lam[1] -= l1g                    # -(1+u) log(1+g)
    rest = [mp.mpf(0)] + Lam[1:]
    E = sexpm1_rest(rest)                            # exp(Lam - Lam0) - 1
    em1 = mp.expm1(Lam[0]); e0 = em1 + 1
    br = [em1 + e0*E[0]] + [e0*E[k] for k in range(1, K+1)]   # (1+g)^(-1-u) S - 1
    pw = [Lg**k/mp.factorial(k) for k in range(K+1)]          # g^u
    tot = smul(br, pw)
    return [v/g for v in tot]

def Mk(alpha, k=K):
    cosa = mp.cos(alpha)
    f = lambda g: integrand_series(g, cosa)[k]
    return mp.quad(f, [0, mp.mpf(1)/8, mp.mpf(1)/2, mp.mpf(7)/8, 1])

if __name__ == '__main__':
    dps = int(sys.argv[1]); mp.mp.dps = dps
    for a in sys.argv[2:]:
        al = mp.mpmathify(eval(a, {'pi': mp.pi, 'mpf': mp.mpf}))
        v = Mk(al)
        print(f"alpha={a}  c={mp.nstr(mp.cos(al/2), 20)}  M_5 = {mp.nstr(v, dps-5)}", flush=True)
