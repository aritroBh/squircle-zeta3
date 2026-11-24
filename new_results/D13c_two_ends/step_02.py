
import json, os
from mpmath import mp, mpf, mpc, pi, factorial, sqrt, cos, sin, exp, log, nstr, zeta, cosh, quad

mp.dps = 80
DIR2 = "../D11b2_stokes_jump"

# ---------- load f_n (450 digits, f_0..f_180) ----------
raw = json.load(open(os.path.join(DIR2,"step2_num450_N180.json")))
print("JSON keys:", list(raw.keys())[:10])
fkey = [k for k in raw.keys() if "f" in k.lower()][0]
print("using key:", fkey, "len:", len(raw[fkey]))
f = [mpf(v) for v in raw[fkey]]
print("f_0..f_6 loaded:")
for n in range(7):
    print("  f_%d = %s" % (n, nstr(f[n], 25)))

# exact zeta forms from D11
fex = {
 0: mpf(1), 1: mpf(-1),
 2: zeta(2)**2/4*0 + zeta(2)/4 + 1,   # zeta2/4 + 1
 3: -zeta(2)/4 - zeta(3)/4 - 1,
 4: -7*zeta(2)**2/32 + zeta(2)/4 + zeta(3)/4 + zeta(4)/32 + 1,
 5: 7*zeta(2)**2/32 + 7*zeta(2)*zeta(3)/16 - zeta(2)/4 - zeta(3)/4 - zeta(4)/32 + 3*zeta(5)/16 - 1,
}
print("\nexact-form check:")
for n,v in fex.items():
    print("  n=%d residual = %s" % (n, nstr(abs(f[n]-v), 4)))

# ---------- load M_k table ----------
def load_mk(fn):
    out = {}
    for line in open(fn):
        p = line.split("\t")
        if len(p) < 2: continue
        k = int(p[0]); mid = p[1].strip().split(" ")[0].strip("[]")
        out[k] = mpf(mid)
    return out
Mk = load_mk("../D13b_Mk_constant/mk_deg9.tsv")
print("\nM_k table: k range", min(Mk), "..", max(Mk), " count", len(Mk))
for k in [0,1,2,10,50,100,400]:
    print("  M_%d = %s" % (k, nstr(Mk[k], 20)))

# ---------- check identity F(eta) = sqrt2/(1+eta) + sqrt2 sum M_k (1+eta)^k ----------
# i.e. M_k = (sqrt2/2)[eps^k]F(eps-1), with index shift.  Compare against quadrature F.
def F_direct(eta, dps_quad=70):
    old = mp.dps; mp.dps = dps_quad
    try:
        g = lambda t: (t*(1-t))**(eta/2) * sqrt(cosh(eta*log(t/(1-t))))
        v = quad(g, [0, mpf('0.25'), mpf('0.5'), mpf('0.75'), 1])
    finally:
        mp.dps = old
    return v

sq2 = sqrt(2)
print("\nidentity check F(eta) vs sqrt2/(1+eta) + sqrt2*sum M_k (1+eta)^k:")
for etas in ["-0.5","-0.3","-0.7","-0.2","-0.9"]:
    eta = mpf(etas)
    s = sq2/(1+eta) + sq2*sum(Mk[k]*(1+eta)**k for k in range(0,401))
    Fd = F_direct(eta)
    print("  eta=%s: series=%s  quad=%s  resid=%s" % (etas, nstr(s,15), nstr(Fd,15), nstr(abs(s-Fd),3)))
