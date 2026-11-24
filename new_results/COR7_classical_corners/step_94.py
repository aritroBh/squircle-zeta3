import sys, json
sys.path.insert(0, '.')
import mpmath as mp

def Li(s, x): return mp.polylog(s, x)
for dps in (300, 420, 600, 900):
    mp.mp.dps = dps
    b = mp.sqrt(2)-1; rt2 = mp.sqrt(2)
    Llog = mp.log(1+rt2); lam = mp.log(2)
    z3, pi2 = mp.zeta(3), mp.pi**2
    W3b = mp.mpf(open("b21_W3_hpl_1150.txt").read())
    Li3b, Li3mb, Li3s = Li(3, b), Li(3, -b), Li(3, 1/rt2)
    Li2b = Li(2, b)
    r3b = json.load(open("b21_W3_relation.json"))
    a03, cf3 = r3b["a0"], r3b["coeffs"]
    B3b = [z3, Li3b, Li3mb, Li3s, Llog*Li2b, lam*Li2b, pi2*Llog, pi2*lam, lam**3, lam**2*Llog, lam*Llog**2, Llog**3]
    nb3 = len(B3b)
    W3bcf = -sum(mp.mpf(cf3[i])*B3b[i] + mp.mpf(cf3[i+nb3])*rt2*B3b[i] for i in range(nb3))/a03
    print(f"dps {dps}: |W3 - W3cf| = {mp.nstr(abs(W3b - W3bcf), 3)}")
