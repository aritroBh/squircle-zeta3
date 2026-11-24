"""verify_closed.py <dps>
Check of the claimed closed forms for c_2, c_1, c_0 (and delta) of W_4(k) = c_0 + c_1 k + c_2 k^2 (k >= 4) at dps higher than the author's 600.
 * numerical c_j, delta: sum of the two half-integrals computed by verify_const.py (my integrand evaluation, exact table from verify_lib)
 * closed forms: transcribed BY HAND from the task text / identify_w4.md (dictionary below), atoms evaluated from scratch with MY OWN
   power series for Li_s(z) (|z| <= 0.7072, cross-checked against mpmath.polylog), constants pi, zeta(3), log 2, L = log(1+sqrt2).
 * also compares my transcription with the author's stored exact coordinates results/id_lll_normalised.json (exact Fractions),
   evaluates delta from that file, and tests c_0 + 4 c_1 + 16 c_2 = W_4(4) against W4pp_800.txt / W4pp_1060.txt.
PRECISION: mp.mp.dps set after all imports and printed with every residual."""
import sys, os, json
sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import mpmath as mp

dps = int(sys.argv[1])
OUT = os.environ["VERIFY_OUT"]
W4 = "<project>/research/08_theory/WTS_weight_k_structure/w4"
mp.mp.dps = dps                                   # after all imports
print("mp.mp.dps =", mp.mp.dps, flush=True)

F = Fr
# ----- hand transcription: label -> (p, q) meaning coefficient p + q*sqrt2 -----
C2 = {"pi4": (F(-3, 960), 0), "pi2 lam^2": (F(-20, 960), 0), "lam^4": (F(-20, 960), 0), "z3 lam": (F(-120, 960), 0)}
C1 = {"lam^4": (F(1, 12), 0), "lam^3 L": (F(-1, 12), 0), "lam^2 L^2": (F(-1, 4), 0), "lam L^3": (F(1, 6), 0), "L^4": (F(1, 3), 0),
      "pi2 lam^2": (F(11, 96), 0), "pi2 lam L": (F(5, 24), 0), "pi2 L^2": (F(-5, 24), 0), "pi4": (F(11, 480), 0),
      "z3 lam": (F(3, 4), 0), "z3 L": (F(-7, 16), 0),
      "Li3b lam": (-1, 0), "Li3b L": (-2, 0), "Li3mb lam": (-1, 0), "Li3mb L": (2, 0), "Li3h L": (4, 0)}
C0 = {"L^4": (F(-17, 8), F(-1, 12)), "lam L^3": (F(-11, 6), F(-7, 3)), "lam^2 L^2": (F(5, 2), F(1, 8)), "lam^3 L": (F(11, 12), F(13, 12)),
      "lam^4": (F(-37, 96), F(13, 48)), "pi2 L^2": (F(89, 48), F(11, 48)), "pi2 lam L": (F(-109, 48), F(7, 6)),
      "pi2 lam^2": (F(13, 48), F(-157, 192)), "Li2b L^2": (5, 1), "Li2b lam L": (-4, 0), "Li2b lam^2": (0, F(-7, 4)),
      "pi4": (F(853, 11520), F(-391, 2880)), "pi2 Li2b": (F(-19, 24), F(-1, 6)), "Li2b^2": (5, 0),
      "z3 L": (F(119, 32), F(3, 8)), "z3 lam": (F(-15, 16), F(77, 64)),
      "Li3b L": (18, 2), "Li3b lam": (4, 6), "Li3mb L": (-18, 2), "Li3mb lam": (4, -6), "Li3h L": (-34, 0), "Li3h lam": (0, -11),
      "Li4(1/2)": (F(-37, 4), F(-1, 8)), "Li4(b)": (1, -104), "Li4(-b)": (-1, -104), "Li4(1/sqrt2)": (0, 2),
      "Li4(-b^2)": (0, F(41, 4)), "Li4(2-sqrt2)": (-10, 16), "Li4(1-1/sqrt2)": (10, 16)}
n_nonzero = sum((1 if F(p) != 0 else 0) + (1 if F(q) != 0 else 0) for p, q in C0.values())
print("c_0 transcription: number of nonzero (p,q) entries =", n_nonzero, "(claimed 51);  c_1:", len(C1), "(claimed 16);  c_2:", len(C2), "(claimed 4)")

# ----- exact comparison with the author's stored coordinates -----
norm = json.load(open(os.path.join(W4, "results", "id_lll_normalised.json")))


def stored(nm):
    return {k: F(v) for k, v in next(iter(norm[nm].values())).items()}


def expand(C):
    out = {}
    for lab, (p, q) in C.items():
        if F(p) != 0: out[lab] = F(p)
        if F(q) != 0: out[lab + " *s2"] = F(q)
    return out


for nm, C in (("c2", C2), ("c1", C1), ("c0", C0)):
    print(f"my hand transcription of {nm} == author's stored exact coordinates (all runs identical):", expand(C) == stored(nm),
          "  all four stored runs identical:", all({k: F(v) for k, v in r.items()} == stored(nm) for r in norm[nm].values()))

# ----- atoms from scratch -----
r2 = mp.sqrt(2); b = r2 - 1; L = mp.log(1 + r2); lam = mp.log(2); pi2 = mp.pi ** 2; z3 = mp.zeta(3)


def Li(s, z):
    z = mp.mpf(z); tot = mp.mpf(0); p = mp.mpf(1); eps = mp.mpf(10) ** (-(mp.mp.dps + 10)); n = 0
    while True:
        n += 1
        p *= z
        term = p / mp.mpf(n) ** s
        tot += term
        if abs(term) < eps * (1 if n < 5 else 1e-0) and n > 5:
            break
    return tot


args4 = {"Li4(1/2)": mp.mpf(1) / 2, "Li4(b)": b, "Li4(-b)": -b, "Li4(1/sqrt2)": 1 / r2, "Li4(-b^2)": -b * b, "Li4(2-sqrt2)": 2 - r2,
         "Li4(1-1/sqrt2)": 1 - 1 / r2}
Li2b = Li(2, b); Li3b = Li(3, b); Li3mb = Li(3, -b); Li3h = Li(3, 1 / r2)
A = {"lam^4": lam ** 4, "lam^3 L": lam ** 3 * L, "lam^2 L^2": lam ** 2 * L ** 2, "lam L^3": lam * L ** 3, "L^4": L ** 4,
     "pi2 lam^2": pi2 * lam ** 2, "pi2 lam L": pi2 * lam * L, "pi2 L^2": pi2 * L ** 2,
     "Li2b lam^2": Li2b * lam ** 2, "Li2b lam L": Li2b * lam * L, "Li2b L^2": Li2b * L ** 2,
     "pi4": pi2 ** 2, "pi2 Li2b": pi2 * Li2b, "Li2b^2": Li2b ** 2, "z3 lam": z3 * lam, "z3 L": z3 * L,
     "Li3b lam": Li3b * lam, "Li3b L": Li3b * L, "Li3mb lam": Li3mb * lam, "Li3mb L": Li3mb * L, "Li3h lam": Li3h * lam, "Li3h L": Li3h * L}
for nm, z in args4.items():
    A[nm] = Li(4, z)
# cross-check my series against mpmath.polylog at this precision
chk = max(abs(Li(4, z) - mp.polylog(4, z)) for z in args4.values())
chk = max(chk, abs(Li2b - mp.polylog(2, b)), abs(Li3b - mp.polylog(3, b)), abs(Li3mb - mp.polylog(3, -b)), abs(Li3h - mp.polylog(3, 1 / r2)))
print("max |my Li_s series - mpmath.polylog| over all atoms:", mp.nstr(chk, 4), " (mp.dps", mp.mp.dps, ")")


def ev(coords):
    tot = mp.mpf(0)
    for lab, c in coords.items():
        base = lab.replace(" *s2", "")
        tot += mp.mpf(c.numerator) / c.denominator * A[base] * (r2 if lab.endswith(" *s2") else 1)
    return tot


def half_sum(which):
    v = mp.mpf(0)
    for h in (0, 1):
        d = json.load(open(os.path.join(OUT, f"const_{which}_h{h}_dps{dps}.json")))
        assert d["mp_dps"] == dps
        v += mp.mpf(d["value"])
    return v


num = {w: half_sum(w) for w in ("c2", "c1", "c0", "delta")}
for w in num:
    print(f"numerical {w} (my integrals) = {mp.nstr(num[w], 60)}")
for nm, C, w in (("c2", C2, "c2"), ("c1", C1, "c1"), ("c0", C0, "c0")):
    val = ev(expand(C))
    print(f"RESIDUAL  closed form {nm} (hand transcription) - my quadrature = {mp.nstr(val - num[w], 5)}   (mp.dps {mp.mp.dps})")
dval = ev(stored("delta"))
print(f"RESIDUAL  delta (author's stored 46-term coordinates) - my quadrature (P(3)-W_4(3)) = {mp.nstr(dval - num['delta'], 5)}   (mp.dps {mp.mp.dps})")
print("   delta =", mp.nstr(num["delta"], 50))

# explicit closed form of c_2 and the exact rewriting
c2cf = -(3 * mp.pi ** 4 + 20 * pi2 * lam ** 2 + 20 * lam ** 4 + 120 * z3 * lam) / 960
print(f"RESIDUAL  explicit c_2 formula - my quadrature = {mp.nstr(c2cf - num['c2'], 5)}   (mp.dps {mp.mp.dps})")
W2 = (pi2 + 6 * lam ** 2) / 24; S = z3 / 8 + pi2 * lam / 24 + lam ** 3 / 12
print(f"RESIDUAL  W_2^2 - log2*S - (7/16) zeta(4) - c_2(quadrature) = {mp.nstr(W2 ** 2 - lam * S - mp.mpf(7) / 16 * mp.zeta(4) - num['c2'], 5)}")

# W_4(4) = c0 + 4 c1 + 16 c2 against the 800/1060-digit file of the earlier M4a pipeline
w44 = num["c0"] + 4 * num["c1"] + 16 * num["c2"]
for fn in ("W4pp_800.txt", "W4pp_1060.txt"):
    p = os.path.join(W4, "..", "..", "M4a_decomposition", fn)
    s = "".join(open(p).read().split()).split("=")[-1]
    nd = len(s.replace("-", "").replace(".", "")) - 1
    ref = mp.mpf(s)
    print(f"W_4(4)=c0+4c1+16c2 (mine) - {fn} ({nd} digits) = {mp.nstr(w44 - ref, 5)}   (mp.dps {mp.mp.dps})")
# W_4(3) from the offset, vs the 48-digit constant printed in the paper (thm:M3 / thm:W4)
w43 = num["c0"] + 3 * num["c1"] + 9 * num["c2"] - num["delta"]
print("W_4(3) = P(3) - delta =", mp.nstr(w43, 60))
print("   minus paper's 48-digit value 0.795785059832240811320359819629726762546527413871 :",
      mp.nstr(w43 - mp.mpf("0.795785059832240811320359819629726762546527413871"), 5))
