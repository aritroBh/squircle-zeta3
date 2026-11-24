
# pipe.py -- pipeline phases for ZQ5 (uses w5lib)
import sys, os, json, pickle, itertools, math, time
from fractions import Fraction
import mpmath as mp
import sympy as sp
import flint
import w5lib as W
from w5lib import QD

OUT = "out"
os.makedirs(OUT, exist_ok=True)

def qd_to_json(x): return [str(x.a), str(x.b)]
def qd_from_json(v, D): return QD(Fraction(v[0]), Fraction(v[1]), D)

def fmt_qd(x):
    d = x.D
    def fs(f): return str(f)
    if x.b == 0: return f"({fs(x.a)})"
    if x.a == 0: return f"({fs(x.b)}*sqrt{d})"
    sgn = "+" if x.b > 0 else "-"
    return f"({fs(x.a)}{sgn}{fs(abs(x.b))}*sqrt{d})"

def cand_dir(D): return os.path.join(OUT, f"d{D}")

# ---------- phase 1: enumeration ----------
def phase_enum(D, B=8):
    C = W.letters_info(D)
    t0 = time.time()
    cands = W.enum_exceptional(C, B=B)
    d = len(C['letters'])
    os.makedirs(cand_dir(D), exist_ok=True)
    data = []
    for c in cands:
        data.append(dict(x=qd_to_json(c['x']), vx=c['vx'], v1m=c['v1m'], h=c['h'],
                         xf=float(c['x'].emb(40))))
    json.dump(dict(D=D, B=B, names=C['names'], disp=C['disp'], S=C['Sdesc'],
                   cands=data), open(os.path.join(cand_dir(D), "candidates.json"), "w"), indent=1)
    print(f"D={D}: {len(cands)} exceptional S-units with |x|<1 (B={B}, {time.time()-t0:.1f}s), {C['Sdesc']}")
    return C, cands

def load_cands(D):
    J = json.load(open(os.path.join(cand_dir(D), "candidates.json")))
    rad = {8:2, 5:5, 12:3, 13:13, 29:29, 21:21, 17:17}[D]
    cands = [dict(x=qd_from_json(c['x'], rad), vx=c['vx'], v1m=c['v1m'], h=c['h']) for c in J['cands']]
    return J, cands

def print_cands(D, cands, C):
    print(f"--- D={D} candidates ({len(cands)}), letters {C['disp']} ---")
    for i, c in enumerate(cands):
        print(f"[{i:2d}] x={fmt_qd(c['x']):>22s} ~ {float(c['x'].emb(30)):+.15f}  vx={c['vx']}  1-x={c['v1m']}  h={c['h']}")

# ---------- phase 2: symbols & kernel ----------
def phase_kernel(D, C, cands):
    d = len(C['letters'])
    N = len(cands)
    A = [W.symbol_li5_dense(c['vx'], c['v1m'], d) for c in cands]
    P, mskeys = W.sym5_basis(d)
    t0 = time.time()
    Wann = W.int_left_kernel([list(col) for col in zip(*P)])  # right kernel of P
    Wm = sp.Matrix(Wann)             # (d^5-m) x d^5
    Am = sp.Matrix(A)                # N x d^5
    Apr = (Am * Wm.T)                # N x (d^5-m)
    K = W.int_left_kernel([list(r) for r in Apr.tolist()])
    # certificates
    Km = sp.Matrix(K)
    assert Km * Apr == sp.zeros(len(K), Apr.shape[1]), "kernel check failed"
    inv = W.snf_invariants(K)
    nullity = Apr.shape[0] - Apr.rank()
    print(f"D={D}: d={d} N={N}  tensor dim={d**5}  Sym^5 dim={len(P)}  quotient dim={len(Wann)}")
    print(f"   rank(symbol map mod Sym^5)={Apr.rank()}  kernel dim={len(K)} (Q-nullity {nullity})  snf={inv}  ({time.time()-t0:.0f}s)")
    json.dump(dict(K=K), open(os.path.join(cand_dir(D), "kernel.json"), "w"))
    return K

# ---------- phase 3: numerics ----------
def phase_numerics(D, C, cands, P=1150, recheck=True):
    d = len(C['letters'])
    t0 = time.time()
    mp.mp.dps = P + 30
    vals = {}
    # logs of letters
    logs = [mp.log(abs(g.emb(P+30))) for g in C['letters']]
    for n, l in zip(C['names'], logs): vals["log:" + n] = l
    vals["pi2"] = mp.pi**2
    vals["z5"] = mp.zeta(5)
    # Li5 values
    li5 = []
    for i, c in enumerate(cands):
        x = c['x'].emb(P+30)
        v = W.li5_smart(x, P)
        li5.append(v)
        vals[f"li5:{i}"] = v
    # cross-check with polylog
    mp.mp.dps = P + 30
    maxdiff = mp.mpf(0)
    for i, c in enumerate(cands):
        x = c['x'].emb(P+30)
        v2 = mp.polylog(5, x)
        maxdiff = max(maxdiff, abs(v2 - li5[i]))
    # L values at two precisions for stability
    L1 = W.Lchi_hurwitz(D, P)
    L2 = W.Lchi_hurwitz(D, P + 150)
    Ldiff = L1 - L2
    vals["L"] = L1
    rad = C['letters'][0].D
    vals["T"] = L1 / mp.sqrt(rad)
    pickle.dump({k: mp.nstr(v, P + 25) for k, v in vals.items()}, open(os.path.join(cand_dir(D), f"vals_{P}.pkl"), "wb"))
    print(f"D={D}: numerics at {P} digits done ({time.time()-t0:.0f}s); "
          f"max|series-polylog|={mp.nstr(maxdiff,3)}; L(P)-L(P+150)={mp.nstr(Ldiff,3)}")
    return vals

def load_vals(D, P=1150):
    raw = pickle.load(open(os.path.join(cand_dir(D), f"vals_{P}.pkl"), "rb"))
    old = mp.mp.dps; mp.mp.dps = P + 100   # parse strings without losing digits
    out = {k: mp.mpf(v) for k, v in raw.items()}
    mp.mp.dps = old
    return out

# ---------- products ----------
def prod_values(C, vals):
    prods = W.weight5_products(C['names'])
    out = []
    for kind, ms, name in prods:
        if kind == "log5":
            v = mp.mpf(1)
            for i in ms: v *= vals["log:" + C['names'][i]]
        elif kind == "pi2log3":
            v = vals["pi2"]
            for i in ms: v *= vals["log:" + C['names'][i]]
        else:
            v = vals["pi2"]**2 * vals["log:" + C['names'][ms[0]]]
        out.append((kind, ms, name, v))
    return out

# ---------- LLL ----------
def lll_rows(vals_list, P):
    mp.mp.dps = P + 30
    n = len(vals_list)
    Cc = mp.power(10, P)
    B = []
    for i, v in enumerate(vals_list):
        row = [0]*(n+1); row[i] = 1; row[n] = int(mp.nint(v*Cc))
        B.append(row)
    M = flint.fmpz_mat(B)
    t0 = time.time()
    R = M.lll(gram='exact')
    rows = [[int(R[i, j]) for j in range(n+1)] for i in range(n)]
    return rows, time.time()-t0

def true_residual(a, vals_list, P):
    mp.mp.dps = 1200   # vals are stored at ~1150 digits; keep full precision
    return sum(mp.mpf(ai)*v for ai, v in zip(a, vals_list))

def phase_lll(D, C, cands, K, vals, P=830, use_atoms=False):
    """vals from load_vals; returns list of relation rows (a_0!=0) and zero rows."""
    mp.mp.dps = 1180   # build combos/products at full stored precision
    prods = prod_values(C, vals)
    N = len(cands)
    li5v = [vals[f"li5:{i}"] for i in range(N)]
    if use_atoms:
        combo_names = [f"Li5({fmt_qd(c['x'])})" for c in cands]
        combo_vals = li5v
        Kmat = None
    else:
        Kmat = K
        combo_names = [f"K{j}" for j in range(len(K))]
        combo_vals = []
        for row in K:
            combo_vals.append(sum(mp.mpf(c)*v for c, v in zip(row, li5v)))
    # drop exact-zero combos
    keep = []
    zero_combos = []
    for nm, v in zip(combo_names, combo_vals):
        if abs(v) < mp.mpf(10)**(-(P-80)):
            zero_combos.append((nm, v))
        else:
            keep.append((nm, v))
    names = ["T", "z5"] + [k[0] for k in keep] + [p[2] for p in prods]
    vals_list = [vals["T"], vals["z5"]] + [k[1] for k in keep] + [p[3] for p in prods]
    rows, dt = lll_rows(vals_list, P)
    rels, zeros = [], []
    for row in rows:
        a = row[:-1]
        mx = max(abs(x) for x in a)
        if mx == 0: continue
        res = true_residual(a, vals_list, P)
        if abs(res) < mp.mpf(10)**(-(P-50)):
            (rels if a[0] != 0 else zeros).append((mx, a, res))
    rels.sort(key=lambda t: t[0])
    zeros.sort(key=lambda t: t[0])
    # EXACT rows: residual near the stored-value precision (spurious LLL rows have
    # residual ~ max|a| * 1e-P, exact ones ~ max|a| * 1e-1150)
    exact_rels  = [t for t in rels  if abs(t[2]) < mp.mpf('1e-1080')]
    exact_zeros = [t for t in zeros if abs(t[2]) < mp.mpf('1e-1080')]
    print(f"D={D} LLL P={P} use_atoms={use_atoms}: dim={len(vals_list)}, {dt:.0f}s; "
          f"relations(a0!=0)={len(rels)}, zero-rels={len(zeros)}, dropped exact-zero combos={len(zero_combos)}")
    for nm, v in zero_combos:
        print(f"   exact-zero combo {nm} ~ {mp.nstr(v,3)}")
    return dict(names=names, vals=vals_list, rels=rels, zeros=zeros,
                exact_rels=exact_rels, exact_zeros=exact_zeros, prods=prods,
                keep=keep, Kmat=Kmat, N=N)
