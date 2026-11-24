
import mpmath as mp, json
mp.mp.dps = 80
h1 = json.load(open('<project>/research/10_frontier/H1_octahedron_end/hi_coeffs_n3_K80.json'))
tau = [mp.mpf(x) for x in h1['eta']]
d   = [mp.mpf(x) for x in h1['d']]
A3 = 2/mp.sqrt(mp.pi**2 + mp.log(2)**2)
th_tau = mp.arg((-mp.log(2)+mp.pi*1j)/2)

# (1) max-envelope root test: E_k = max_{j in [k-3,k+3]} (|a_j|/j!)^(1/j)
def envelope_root(a, rng):
    out = []
    for k in rng:
        lo = max(2, k-3); hi = min(len(a)-1, k+3)
        out.append((k, max((abs(a[j])/mp.factorial(j))**(1/j) for j in range(lo,hi+1))))
    return out
print("max-envelope root test (|a_k|/k!)^(1/k), predicted ->", mp.nstr(A3,8))
for k, v in envelope_root(tau, [20,30,40,50,60,70,76]):
    print(f"  tau k={k}: {mp.nstr(v,9)}")
print("d-series (predicted -> 1/|d*| = 1.286216):")
for k, v in envelope_root(d, [20,30,40,50,60,70,76]):
    print(f"  d   k={k}: {mp.nstr(v,9)}")

# (2) model fit: a_k/(k!) = k^beta (C1 cos(k th) + C2 sin(k th)) with th known; grid over A
def fit_A(a, th, klo, khi):
    best = None
    for i in range(-40, 41):
        A = A3*(1 + i*mp.mpf('0.002'))
        ks = list(range(klo, khi+1))
        # for each beta in small grid solve linear C1,C2 and get residual
        for j in range(-30, 11):
            beta = j*mp.mpf('0.1')
            S = [[mp.mpf(0)]*2 for _ in range(2)]; T=[mp.mpf(0)]*2
            for k in ks:
                x1 = k**beta*mp.cos(k*th); x2 = k**beta*mp.sin(k*th); y = a[k]/mp.factorial(k)/A**k
                S[0][0]+=x1*x1; S[0][1]+=x1*x2; S[1][1]+=x2*x2; T[0]+=x1*y; T[1]+=x2*y
            S[1][0]=S[0][1]
            det = S[0][0]*S[1][1]-S[0][1]**2
            if abs(det) < mp.mpf('1e-60'): continue
            C1 = (T[0]*S[1][1]-T[1]*S[0][1])/det; C2=(S[0][0]*T[1]-S[0][1]*T[0])/det
            res = mp.mpf(0)
            for k in ks:
                x1=k**beta*mp.cos(k*th); x2=k**beta*mp.sin(k*th)
                r = a[k]/mp.factorial(k)/A**k - C1*x1 - C2*x2
                res += r*r
            res = mp.sqrt(res/len(ks))
            scale = mp.sqrt(C1*C1+C2*C2)
            if best is None or res/scale < best[0]:
                best = (res/scale, A, beta, C1, C2)
    return best
relres, Afit, betafit, C1, C2 = fit_A(tau, th_tau, 40, 79)
print(f"\n(2) grid fit tau: A={mp.nstr(Afit,9)} (pred {mp.nstr(A3,9)}), beta={mp.nstr(betafit,5)}, rel resid={mp.nstr(relres,3)}, |Z|={mp.nstr(mp.sqrt(C1**2+C2**2),6)}")
