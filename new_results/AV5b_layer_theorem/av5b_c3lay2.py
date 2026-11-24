
import mpmath as mp
import json
mp.mp.dps = 70

def load_prof(mm):
    return [(mp.mpf(s), mp.mpf(a), mp.mpf(b), mp.mpf(c)) for s,a,b,c in json.load(open(f"ckpt_prof_m{mm}.json"))]
def trapz(ss, ds, vals):
    return mp.fsum((vals[i]+vals[i+1])/2*ds for i in range(len(ss)-1))

for mm in [2, 3]:
    prof = load_prof(mm)
    ss = [p[0] for p in prof]; ds = ss[1]-ss[0]
    H = mp.sin(mp.pi/(2*mm)); alpha = mp.tan(mp.pi/(2*mm))
    a1 = mp.mpf(0) if mm==2 else mp.log(2)
    g1 = [p[1]-a1 for p in prof]; g2 = [p[2] for p in prof]
    I_g2   = trapz(ss, ds, g2)
    I_g1sq = trapz(ss, ds, [x**2/2 for x in g1])
    I_sg1  = trapz(ss, ds, [(alpha*abs(ss[i])+a1)*g1[i] for i in range(len(ss))])
    c3 = H*(I_g2 + I_g1sq + I_sg1)
    print(f"m={mm}: c3_lay = {mp.nstr(c3, 28)}")
    if mm == 2:
        print("   target sqrt2 z3/4:", mp.nstr(mp.sqrt(2)*mp.zeta(3)/4, 28))
        print("   pieces: I_g1sq =", mp.nstr(I_g1sq, 20), " (z3/8 =", mp.nstr(mp.zeta(3)/8, 20), ")  I_sg1 =", mp.nstr(I_sg1, 20), " (3z3/8 =", mp.nstr(3*mp.zeta(3)/8, 20), ")")
    if mm == 3:
        z3 = mp.zeta(3); ln2 = mp.log(2); pi2 = mp.pi**2
        basis = [z3, z3*ln2, pi2*ln2, ln2**3, pi2, ln2**2, ln2, mp.mpf(1)]
        M = [c3] + basis
        rel = mp.pslq(M, tol=mp.mpf('1e-18'), maxcoeff=10**8, maxsteps=10000)
        print("   PSLQ:", rel)
        if rel: print("    resid:", mp.nstr(mp.fsum(mp.mpf(r)*x for r,x in zip(rel,M)),6))
        print("   target C3 z3/36:", mp.nstr(-9*mp.sqrt(3)/8*z3/36, 28))
        print("   pieces: I_g2 =", mp.nstr(I_g2, 18), " I_g1sq =", mp.nstr(I_g1sq, 18), " I_sg1 =", mp.nstr(I_sg1, 18))
