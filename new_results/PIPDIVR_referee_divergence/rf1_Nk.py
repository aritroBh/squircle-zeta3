# rf1: N_k (k = 0..K) from the s-variable moment expansion, in arb ball arithmetic.
# pi_p/4 = 2^{-h} + Ehat,  Ehat = h*kap * int_0^inf G(kap*s,h) [(1+e^{-s})^h - 1] ds,
# G(y,h) = e^{-hy}(1+e^{-y})^{-1-h} = sum_j g_j(h) y^j,  (1+e^{-s})^h - 1 = sum_{r>=1} h^r L^r / r!,
# so termwise  Ehat ~ sum_j h^{j+2}(1-h)^{-j-1} g_j(h) sum_r h^r M_{j,r}/r!,  M_{j,r} = int s^j L^r ds.
# Usage: rf1_Nk.py K PREC S OUTFILE
import sys, time
from flint import arb, acb, ctx
from rf_series import mul, log, exp, inv

K = int(sys.argv[1]); PREC = int(sys.argv[2]); S = int(sys.argv[3]); OUT = sys.argv[4]
ctx.prec = PREC
t0 = time.time()
N = K + 2  # series length

# --- moments M_{j,r}, j + r <= K - 1, rigorous balls (integral on [0,S] + tail bound) ---
def make_f(j, r):
    def f(s, analytic):
        L = (1 + (-s).exp()).log(analytic=analytic)
        return s**j * L**r
    return f
M = {}
rel = arb(2) ** (-(PREC - 8)); ab = arb(2) ** (-(PREC + 2000))
for j in range(0, K):
    for r in range(1, K - j):
        assert S >= 2 * j
        v = acb.integral(make_f(j, r), 0, S, rel_tol=rel, abs_tol=ab).real
        # tail: int_S^inf s^j L^r <= int_S^inf s^j e^{-rs} <= 2 S^j e^{-rS}/r  (S >= 2j)
        tail = 2 * arb(S) ** j * (-arb(r * S)).exp() / r
        v = v + arb(0, tail.upper())
        M[(j, r)] = v
print("moments done", len(M), "%.1fs" % (time.time() - t0), flush=True)

# --- bivariate coefficients g[j][i] = [y^j h^i] G(y,h) ---
fac = [arb(1)]
for n in range(1, 2 * N + 2):
    fac.append(fac[-1] * n)
emy = [arb(-1) ** n / fac[n] for n in range(N)]          # e^{-y}
u = [emy[0] + 1] + emy[1:]                                  # 1 + e^{-y}
q = log(u, N, arb(2).log())                                 # log(1+e^{-y})
P = inv(u, N)                                               # e^{-q}
B = [-q[0], -q[1] - 1] + [-c for c in q[2:]]                # -(q + y)
g = [[None] * N for _ in range(N)]                          # g[j][i]
Bi = [arb(1)] + [arb(0)] * (N - 1)
for i in range(N):
    term = mul(P, Bi, N)
    for j in range(N):
        g[j][i] = term[j] / fac[i]
    Bi = mul(Bi, B, N)

# --- assemble [h^{k+1}] Ehat ---
Ecoef = [arb(0)] * (K + 2)
def binom(a, b):
    return fac[a] / (fac[b] * fac[a - b])
for j in range(0, K - 1):
    top = K + 1 - (j + 2)              # need h-powers up to this inside the bracket
    if top < 1:
        continue
    om = [binom(j + i, i) for i in range(top + 1)]      # (1-h)^{-j-1}
    gj = [g[j][i] for i in range(top + 1)]
    Jj = [arb(0)] + [M[(j, r)] / fac[r] if (j, r) in M else arb(0) for r in range(1, top + 1)]
    prod = mul(mul(om, gj, top + 1), Jj, top + 1)
    for e in range(top + 1):
        Ecoef[j + 2 + e] += prod[e]
l2 = arb(2).log()
with open(OUT, "w") as fh:
    for k in range(0, K):
        Nk = (-l2) ** (k + 1) / fac[k + 1] + Ecoef[k + 1]
        fh.write("%d %s\n" % (k, Nk.str(radius=True, more=True, condense=0) if False else Nk.str(60, radius=True)))
print("done %.1fs" % (time.time() - t0))
