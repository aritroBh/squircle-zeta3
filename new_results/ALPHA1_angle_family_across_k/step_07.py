"""step_07.py -- SQUIRCLE CONTROL: c = 1/sqrt2 (r = 3-2sqrt2, QSI field), k = 1..8.
My w-kernels, mapped by v = a w, must equal the squircle kern_k.pkl (exact).
Real residues: res_{w=rho/a} = res_{v=rho}.  theta-parts: coeff_j(w) = a^{1-j} coeff_j(v).
"""
import sys, time, pickle
sys.path.insert(0, '.')
from fractions import Fraction as Fr
from cwfield import QSI
from cwintegrand import CornerIntegrand
from cwred import Reducer
from cwbuild import build_integrand_pf

TH = '<project>/research/08_theory/WTS_weight_k_structure/theory/'
r_frac = None
c = QSI(0, Fr(1,2))          # 1/sqrt2 = sqrt2/2
r = QSI(3, -2)               # 3 - 2 sqrt2 = b^2
a = QSI(1, 1)                # 1 + sqrt2
F = QSI
red = Reducer(F, r)
ci = CornerIntegrand(F, c, r)

def qsi_from_plain(t):
    return QSI(*[Fr(x) for x in t])

def qsi2(t):   # (p0,p1) Fractions -> QSI (real)
    return QSI(Fr(t[0]), Fr(t[1]), 0, 0)

def load_squircle(k):
    kd = pickle.load(open(TH + f'kern_{k}.pkl', 'rb'))
    out = {}
    for key, P in kd.items():
        real = {qsi2(rho): qsi_from_plain(cf) for rho, cf in P['real'].items()}
        ipos = [qsi_from_plain(e) for e in P['ipos']]
        ineg = [qsi_from_plain(e) for e in P['ineg']]
        const = qsi_from_plain(P['const'])
        out[key] = [real, ipos, ineg, const]
    return out

def plain_rho(rho):
    return qsi2(rho)

KMAX = 8
allok = True
for k in range(1, KMAX+1):
    t0 = time.time()
    ipf = build_integrand_pf(ci, red, k)
    fin = red.reduce_k(ipf)
    sq = load_squircle(k)
    ok = True
    msgs = []
    # compare keys present in both (squircle kern has keys (i,j) i+j<=k? check)
    for key, Psq in sq.items():
        if key not in fin:
            ok = False; msgs.append(f"missing key {key}"); continue
        Pw = fin[key]
        # real poles
        for rho_v, cf in Psq[0].items():
            rho_w = rho_v / a
            got = Pw[0].get(rho_w, QSI(0))
            if got != cf:
                ok = False; msgs.append(f"key{key} real {rho_v}: sq={cf} mine={got}")
        for rho_w, cf in Pw[0].items():
            rho_v = rho_w * a
            if rho_v not in Psq[0] and not cf.is_zero():
                ok = False; msgs.append(f"key{key} extra real pole {rho_w}")
        # theta parts
        ipos_sq = Psq[1]
        ipos_mine_scaled = [ (a**(1-j))*e for j, e in enumerate(Pw[1], start=1) ]
        if len(ipos_sq) != len(ipos_mine_scaled) or any(x != y for x, y in zip(ipos_sq, ipos_mine_scaled)):
            ok = False; msgs.append(f"key{key} theta+: sq={ipos_sq} mine_scaled={ipos_mine_scaled}")
    nkeys = len(sq)
    print(f"k={k}: squircle keys={nkeys}, match={'OK' if ok else 'FAIL'} ({time.time()-t0:.1f}s)")
    for m in msgs[:6]: print("   ", m)
    allok = allok and ok
print("SQUIRCLE CONTROL:", "PASS" if allok else "FAIL")
