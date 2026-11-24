# 00_reference.py -- independent 100-digit quadrature values of P1, P2a, P2c (integrals over s in [0,1])
import mpmath as mp, json
mp.mp.dps = 110

def qs(s): return mp.sqrt(s*s + (1 - s)**2)

def p1(s):
    m = s*(1 - s)
    extra = int(max(0, -mp.log10(m))) + 10 if m > 0 else 0
    with mp.workdps(mp.mp.dps + extra):
        q = mp.sqrt(1 - 2*m)
        v = (2*mp.log(1 + q) - (1 + q)*mp.log(2*q))/(2*m)     # 1-q^2 = 2s(1-s)
    return +v
def p2a(s):
    q = qs(s); return mp.asinh(q)/mp.sqrt(1 + q*q)
def p2c(s):
    q = qs(s); return mp.asinh(1/q)/mp.sqrt(1 + q*q)

out = {}
for name, f in (("P1", p1), ("P2a", p2a), ("P2c", p2c)):
    v = 2*mp.quad(f, [0, mp.mpf(1)/8, mp.mpf(1)/4, mp.mpf(1)/2])
    out[name] = mp.nstr(v, 100)
    print(name, "=", mp.nstr(v, 100))
Q0 = 6*(mp.mpf(out["P1"]) - mp.mpf(out["P2a"]) - mp.mpf(out["P2c"]))
out["Q0_from_pieces"] = mp.nstr(Q0, 100)
print("6(P1-P2a-P2c) =", mp.nstr(Q0, 100))
json.dump(out, open("reference_100.json", "w"), indent=1)
