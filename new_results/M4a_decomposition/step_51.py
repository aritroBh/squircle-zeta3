import sys, os, pickle
sys.dont_write_bytecode = True
sys.path.insert(0, ".")
sys.path.insert(0, "<project>/research/08_theory/A3_W4_proof")
import mpmath as mp
from flint import acb, ctx
from w4_symbol import conjecture_terms
import dilog_exact as DE

mp.mp.dps = 1250
ctx.prec = 4200
BASE = "<project>"
M3file = acb(open(os.path.join(BASE, "paper/squircle_note/outputs/M3_methodA_1200.txt")).read().strip())

tails = pickle.load(open("m3_tails_exact.pkl", "rb"))
base_fl = {"lam": acb(2).log(), "L": (acb(2).sqrt()+1).log(), "pi2": acb.pi()**2, "z3": acb(3).zeta()}
def pval_fl(P):
    tot = acb(0)
    for k, c in P.d.items():
        x = acb(c.p.numerator)/c.p.denominator + acb(2).sqrt()*acb(c.q.numerator)/c.q.denominator
        for a in k:
            if isinstance(a, str):
                x *= base_fl[a]
            else:
                z = DE.float_exact(a[1])
                x *= acb(mp.nstr(z, 1200)).polylog(acb(2 if a[0]=="Li2" else 3))
        tot += x
    return tot
W3p_fl = pval_fl(tails["W3p"])
T_fl = pval_fl(tails["T"])
W4th_fl = acb(0)
for c, atoms in conjecture_terms():
    t = acb(c.p.numerator)/c.p.denominator + acb(2).sqrt()*acb(c.q.numerator)/c.q.denominator
    for a in atoms:
        if a[0]=="li":
            z = acb(c2p := None) if False else None
    # do it properly below
def atom_fl(a):
    if a[0]=="li":
        zz = a[2]
        zf = mp.mpf(zz.p.numerator)/zz.p.denominator + mp.sqrt(2)*mp.mpf(zz.q.numerator)/zz.q.denominator
        return acb(mp.nstr(zf, 1200)).polylog(acb(a[1]))
    if a[0]=="log":
        zz = a[1]
        zf = mp.mpf(zz.p.numerator)/zz.p.denominator + mp.sqrt(2)*mp.mpf(zz.q.numerator)/zz.q.denominator
        return acb(mp.nstr(zf,1200)).log()
    if a[0]=="pi2": return acb.pi()**2
    if a[0]=="z3": return acb(3).zeta()
W4th_fl = acb(0)
for c, atoms in conjecture_terms():
    t = acb(c.p.numerator)/c.p.denominator + acb(2).sqrt()*acb(c.q.numerator)/c.q.denominator
    for a in atoms:
        t *= atom_fl(a)
    W4th_fl += t

D = M3file - W3p_fl - T_fl - W4th_fl
print("D = M3file - W3' - T - W4th =", D.real.mid().str(30, radius=False, more=True), " rad", D.real.rad().str(5))
print("W4th_fl   =", W4th_fl.real.mid().str(80, radius=False, more=True))
print("M3file-W3'-T =", (M3file-W3p_fl-T_fl).real.mid().str(80, radius=False, more=True))
