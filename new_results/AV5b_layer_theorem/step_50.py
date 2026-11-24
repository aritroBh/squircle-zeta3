import os
work = "<project>/research/08_theory/AV5b_layer_theorem"
os.chdir(work)
step = r'''
import mpmath as mp
import json
mp.mp.dps = 70

def load_prof(mm):
    data = json.load(open(f"ckpt_prof_m{mm}.json"))
    return [(mp.mpf(s), mp.mpf(a), mp.mpf(b), mp.mpf(c)) for s,a,b,c in data]

def c3_lay(mm, a1, a2):
    prof = load_prof(mm)
    H = mp.sin(mp.pi/(2*mm)); alpha = mp.tan(mp.pi/(2*mm))
    # trapezoid over the grid
    ss = [p[0] for p in prof]
    ds = ss[1]-ss[0]
    def integrand(p):
        s, D1, D2, D3 = p
        g1 = D1 - a1; g2 = D2 - a2; g3 = D3
        G2 = g2 + g1**2/2
        G3 = g3 + g1*g2 + g1**3/6
        return G3 + (alpha*abs(s)+a1)*(G2 + a1*g1) + (a2 + a1**2/2)*g1 - (s**2/2)*g1
    vals = [integrand(p) for p in prof]
    I = mp.fsum((vals[i]+vals[i+1])/2*ds for i in range(len(vals)-1))
    return H*I, prof, integrand

for mm, tgt in [(2, mp.sqrt(2)*mp.zeta(3)/4), (3, None)]:
    a1 = mp.mpf(0) if mm==2 else mp.log(2)
    a2 = mp.mpf(0)
    v, prof, integrand = c3_lay(mm, a1, a2)
    print(f"m={mm}: c3_lay = {mp.nstr(v, 30)}")
    if tgt: print(f"   target sqrt2 z3/4 = {mp.nstr(tgt,30)}  resid {mp.nstr(v-tgt,6)}")
    if mm == 3:
        z3 = mp.zeta(3); ln2 = mp.log(2); pi2 = mp.mp.pi**2 if False else mp.pi**2
        basis = [z3, z3*ln2, pi2*ln2, ln2**3, pi2, ln2**2, ln2, mp.mpf(1)]
        M = [v] + basis
        rel = mp.pslq(M, tol=mp.mpf('1e-20'), maxcoeff=10**8, maxsteps=10000)
        print("   PSLQ c3_lay:", rel)
        if rel: print("    resid:", mp.nstr(mp.fsum(mp.mpf(r)*x for r,x in zip(rel,M)),6))
        print("   C3*z3/(12*3) =", mp.nstr(-9*mp.sqrt(3)/8*z3/36, 30))
'''
with open("av5b_c3lay.py","w") as f: f.write(step)
exec(step)
