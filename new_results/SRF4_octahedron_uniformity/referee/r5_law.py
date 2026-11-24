# Referee item 5: validate the coefficients against quadrature, then compare with the law
#   a_k/k! = -sqrt(3 pi/8) k^{-3/2} Re z0^{1/2-k} + o(k^{-3/2} R^{-k}).
import sys, json
import mpmath as mp
from r1_representation import S_T
mp.mp.dps = 40
fn = sys.argv[1]
d = json.load(open(fn)); K = d["K"]
a = [mp.mpf(x) for x in d["a"]]
if len(sys.argv) > 2:  # precision check
    b = [mp.mpf(x) for x in json.load(open(sys.argv[2]))["a"]]
    print("max rel diff between precisions, k<=K:", mp.nstr(max(abs(a[k]-b[k])/abs(a[k]) for k in range(1, K+1)), 3))
z0 = mp.mpc(-mp.log(2), mp.pi)/2; R = abs(z0); th = mp.arg(z0); kap = mp.sqrt(3*mp.pi/8)

# validation: truncated series vs quadrature at small tau (truncation error ~ a_{K+1} tau^{K+1})
if K <= 60:
    mp.mp.dps = 25
    for tau in ['0.02', '0.03']:
        tau = mp.mpf(tau)
        Sq, _ = S_T(tau)
        ser = sum(a[k]*tau**k for k in range(K+1))
        print(f"tau={tau}: S/S(1) quadrature={mp.nstr(Sq/(4*mp.sqrt(3)),22)}  series(k<={K})={mp.nstr(ser,22)}  diff={mp.nstr(Sq/(4*mp.sqrt(3))-ser,3)}  next term ~{mp.nstr(a[K]*tau**K,3)}")
    mp.mp.dps = 40

print("\n k   |a_k/k!|^(1/k)   R|..|^(1/k)  ratio a/law     residual=(a/k!-law)/(k^-1.5 R^-k)   k*residual   |cos((k-1/2)th)|")
rows = []
for k in range(5, K+1):
    al = a[k]/mp.factorial(k)
    law = -kap*mp.power(k, -1.5)*mp.re(z0**(mp.mpf(1)/2-k))
    sc = mp.power(k, -1.5)*R**(-k)
    res = (al-law)/sc
    rows.append((k, al, law, res))
    if k % 5 == 0 or k in (33, 91):
        print(f"{k:3d}  {mp.nstr(abs(al)**(mp.mpf(1)/k),8):>10}  {mp.nstr(R*abs(al)**(mp.mpf(1)/k),6):>8}  {mp.nstr(al/law,6):>10}  {mp.nstr(res,6):>12}  {mp.nstr(k*res,5):>9}  {mp.nstr(abs(mp.cos((k-0.5)*th)),4)}")
print("1/R =", mp.nstr(1/R, 16))

# amplitude/phase fit: a_k/(k! k^{-3/2} R^{-k}) = A cos(k th) + B sin(k th) + sum_{j=1}^{J} (C_j cos + D_j sin)/k^j
def fit(k1, k2, J):
    ks = list(range(k1, k2+1))
    cols = lambda k: [mp.cos(k*th), mp.sin(k*th)] + sum([[mp.cos(k*th)/mp.mpf(k)**j, mp.sin(k*th)/mp.mpf(k)**j] for j in range(1, J+1)], [])
    M = mp.matrix([cols(k) for k in ks])
    y = mp.matrix([a[k]/mp.factorial(k)/(mp.power(k, -1.5)*R**(-k)) for k in ks])
    x = mp.lu_solve(M.T*M, M.T*y)
    A, B = x[0], x[1]
    amp = mp.sqrt(A**2+B**2); ph = mp.atan2(-B, A)  # amp cos(k th + ph)
    return amp, ph
print("\ntheory: amp = sqrt(3 pi R/8) =", mp.nstr(mp.sqrt(3*mp.pi*R/8), 10), "  phase = pi - th/2 =", mp.nstr(mp.pi-th/2, 10))
for (k1, J) in [(15, 2), (15, 3), (20, 3), (12, 4)]:
    amp, ph = fit(k1, K, J)
    print(f"fit k={k1}..{K}, {J} 1/k-orders: amp={mp.nstr(amp,10)}  phase={mp.nstr(ph % (2*mp.pi),10)}")

# two-consecutive floor, each member on its own scale k^{-3/2} R^{-k}
sc = lambda r: abs(r[1])*mp.power(r[0], 1.5)*R**r[0]
for k0 in (10, 20, 30, 60):
    sel = [i for i in range(len(rows)-1) if rows[i][0] >= k0]
    if sel:
        print(f"min over k>={k0} of max(own-scale |a_k|, |a_(k+1)|) =", mp.nstr(min(max(sc(rows[i]), sc(rows[i+1])) for i in sel), 6),
              " (asymptotic floor sin((pi-th)/2)*sqrt(3piR/8) =", mp.nstr(mp.sin((mp.pi-th)/2)*mp.sqrt(3*mp.pi*R/8), 6), ")")

# d-series: b_n = sum_k a_k (-1)^{n-k} C(n-1,k-1); law b_n/n! = -sqrt(3pi/4) n^{-3/2} Im z0^{1/2-n}
print("\n n   |b_n/n!|^(1/n)   residual_d   n*residual_d")
for n in range(5, K+1):
    bn = sum(a[k]*(-1)**(n-k)*mp.binomial(n-1, k-1) for k in range(1, n+1))
    bl = bn/mp.factorial(n)
    law = -mp.sqrt(3*mp.pi/4)*mp.power(n, -1.5)*mp.im(z0**(mp.mpf(1)/2-n))
    res = (bl-law)/(mp.power(n, -1.5)*R**(-n))
    if n % 5 == 0:
        print(f"{n:3d}  {mp.nstr(abs(bl)**(mp.mpf(1)/n),8):>10}  {mp.nstr(res,6):>10}  {mp.nstr(n*res,5):>8}")
