# rf3: N_k are the asymptotic coefficients of pi_p/4: compare optimally truncated series with pi_p computed
# directly from the definition pi_p = (2/p) int_0^1 (t^{1-p}+(1-t)^{1-p})^{1/p} dt (own code).
# With t = v^p on [0,1/2]:  pi_p/4 = 2^{-h} + int_0^{2^{-h}} [(1 + r^{p-1})^h - 1] dv,  r = v^p/(1-v^p).
import sys, re
import mpmath as mp
NF = sys.argv[1]; ps = [int(x) for x in sys.argv[2].split(",")]; mp.mp.dps = int(sys.argv[3])
N = {}
for line in open(NF):
    k, rest = line.split(" ", 1)
    m = re.match(r"\[(\S+) \+/- (\S+)\]", rest.strip())
    N[int(k)] = mp.mpf(m.group(1))
K = max(N) + 1
def pip4(p):
    h = mp.mpf(1)/p; vmax = mp.mpf(2)**(-h)
    def f(v):
        vp = v**p; r = vp/(1-vp)
        return mp.expm1(h*mp.log1p(r**(p-1)))
    a = vmax*mp.exp(-300*h*h)
    pts = [a + (vmax-a)*mp.mpf(j)/40 for j in range(41)]
    return vmax + mp.quad(f, pts), mp.quad(f, [0, a])   # second = the omitted [0,a] piece (must be negligible)
def Ehat_split(p):   # the D12 split, as a cross-check of the direct value
    h = mp.mpf(1)/p; kap = h/(1-h)
    f = lambda s: mp.exp(-h*kap*s)*(1+mp.exp(-kap*s))**(-1-h)*mp.expm1(h*mp.log1p(mp.exp(-s)))
    return mp.mpf(2)**(-h) + h*kap*mp.quad(f, [0, 1, 5, 20, 60, 150, 400, mp.inf])
print("p  |direct - split|  best K  |pi_p/4 - S_K| at best K   e^{-pi p}   ratio   smallest |N_k h^{k+1}| (k)")
for p in ps:
    h = mp.mpf(1)/p
    d, d2 = pip4(p); sp = Ehat_split(p)
    S = mp.mpf(1); best = (mp.inf, None); terms = []
    for k in range(K):
        S += N[k]*h**(k+1); err = abs(d - S); terms.append(abs(N[k]*h**(k+1)))
        if err < best[0]: best = (err, k+1)
    kmin = min(range(K), key=lambda k: terms[k] if k > 3 else mp.inf)
    print(p, mp.nstr(abs(d-sp), 3), best[1], mp.nstr(best[0], 4), mp.nstr(mp.exp(-mp.pi*p), 4),
          mp.nstr(best[0]/mp.exp(-mp.pi*p), 4), mp.nstr(terms[kmin], 4), kmin, "  [0,a] piece:", mp.nstr(d2, 3))
