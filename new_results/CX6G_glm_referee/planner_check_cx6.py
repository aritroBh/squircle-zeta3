# Planner check of CX6 (decay-law rest bound): constants, assembly chain, (13), (14), and claim (2) on exact M_k data.
import mpmath as mp, re, sys
mp.mp.dps = 50
pi = mp.pi
a = pi/(2*mp.sqrt(2)); x = 2*mp.e**(pi/4)
K0 = 64/(pi*a)*mp.e**(a/2+pi/4); T0 = 32/(3*pi*a)*mp.e**(3*a/2+3*pi/4)
Cs = 2 + T0 + K0*(mp.e**(pi/4) + mp.e**x + mp.e**(mp.e*x))
print("K0 =", mp.nstr(K0, 25)); print("T0 =", mp.nstr(T0, 25)); print("C* =", mp.nstr(Cs, 25), " < 1.06e7:", Cs < 10600000)
d = 2 - mp.sqrt(pi); conv = (3/(mp.e*d))**3; tot = mp.mpf('254.23') + 10600000*conv
print("conv =", mp.nstr(conv, 20), " total =", mp.nstr(tot, 20), " < 1.21e9:", tot < 1210000000)
# assembly chain: sum of the explicit pieces (15),(17),(18),J_k versus C* e^{-2 sqrt N}, N = 100..3000
worst = 0
for N in list(range(100, 400)) + list(range(400, 3001, 50)):
    k = N - 1
    s15 = mp.e**(pi/4)*mp.e**(-mp.sqrt(3*pi/2*k))
    s15 += mp.fsum(x**i/mp.factorial(i)*mp.e**(-mp.sqrt(3*pi/4*(i+1)*(k-i))) for i in range(1, k+1))
    tot_pieces = K0*s15 + T0*mp.e**(-mp.sqrt(9*pi/4*k)) + mp.log(2)**N/mp.factorial(N) + 2/(pi*k)*mp.mpf(3)**(-k)
    r = tot_pieces/(Cs*mp.e**(-2*mp.sqrt(N)))
    worst = max(worst, r)
print("assembly: max over N of pieces/(C* e^{-2sqrtN}) =", mp.nstr(worst, 10), "(must be <= 1)")
# (13): min over r of |1+exp(i pi s)| on Gamma, s = 1/2 + r e^{i pi/4}
w = mp.e**(1j*pi/4)
m = min(abs(1+mp.e**(1j*pi*(mp.mpf(1)/2 + r*w))) for r in [mp.mpf(t)/400 for t in range(0, 4001)])
print("(13): min |1+e^{i pi s}| on Gamma, r in [0,10] =", mp.nstr(m, 12), "(claim >= 1/2)")
# (14): F(j,q;lam) by quadrature versus the bound
mp.mp.dps = 30
worst14 = 0
for (j, q, lam) in [(0,1,1),(1,1,1),(5,1,2.2),(20,1,3.3),(50,2,1.1),(100,1,3.33),(200,3,0.7),(400,1,3.33)]:
    q = mp.mpf(q); lam = mp.mpf(lam)
    f = lambda r: mp.e**(-lam*r)*abs(mp.mpf(1)/2 + r*w)**j/abs(mp.mpf(1)/2 + r*w + q)**(j+1)
    F = mp.quad(f, [0, 1, 5, 20, 80, mp.inf])
    B = 4/lam*mp.e**(lam/2 + lam*q/mp.sqrt(2))*mp.e**(-mp.sqrt(3*lam*j*q/mp.sqrt(2)))
    worst14 = max(worst14, F/B)
    print(f"(14) j={j} q={q} lam={lam}: F={mp.nstr(F,8)} bound={mp.nstr(B,8)} ratio={mp.nstr(F/B,6)}")
print("(14): worst ratio =", mp.nstr(worst14, 8), "(must be <= 1)")
# claim (2) on the exact table M_k, k = 0..400 (N = k+1 >= 100)
mp.mp.dps = 90
rows = open(sys.argv[1]).read().split("\n")
mx = 0; kmx = None
for line in rows:
    if not line.strip(): continue
    k, val = line.split("\t")[:2]
    k = int(k); N = k + 1
    if N < 100: continue
    Mk = mp.mpf(re.match(r"\[(\S+)", val).group(1))
    lhs = abs(Mk + mp.sqrt(pi)/N*mp.e**(-mp.sqrt(pi*N))*mp.cos(mp.sqrt(pi*N)))
    ratio = lhs*mp.mpf(N)**1.5*mp.e**(mp.sqrt(pi*N))
    if ratio > mx: mx, kmx = ratio, k
print("claim (2) on table: max |M_k + main| N^{3/2} e^{sqrt(pi N)} =", mp.nstr(mx, 10), "at k =", kmx, "(paper constant 1.21e9; main-term thm 254.23)")
