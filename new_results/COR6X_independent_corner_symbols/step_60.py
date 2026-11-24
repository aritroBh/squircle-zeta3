t = r'''
# step_06e.py -- DECISIVE: weight-4 part of M3(2pi/3) vs the paper's Li4 combination
import sys; sys.path.insert(0, '.')
import mpmath as mp, json, re
from hpl_num import hpl
from lll_rel import find_relations
mp.mp.dps = 200
def load(f):
    data = json.load(open("<project>/research/10_frontier/G1_zeta3_any_angle/" + f))
    out = []
    for s, c in data:
        m = re.match(r"c\^(-?\d+)\s+H\[([^\]]*)\]", s)
        if m:
            out.append((int(m.group(1)), tuple(int(t) for t in m.group(2).split(',')), mp.mpf(c)))
        else:
            out.append((None, None, mp.mpf(c)))
    return out
R3 = load("hpl_fit_k3_(0,0,0).json")
PSI = load("hpl_fit_k3_(1,0,0).json")
def Omega(c0):
    return mp.fsum(cf * (c0 ** e) * hpl(w, c0) for e, w, cf in R3 if w is not None and len(w) == 4)
def Psi_wt2(c0):
    return mp.fsum(cf * (c0 ** e) * hpl(w, c0) for e, w, cf in PSI if w is not None and len(w) == 2)
def Phi_wt1(c):
    return -(1 - c**2)/(2*c) * mp.atanh(c)
c = mp.mpf(1)/2
W4M = Omega(c) - mp.zeta(3)*Phi_wt1(c) + mp.pi**2/6*Psi_wt2(c)
print("wt4(M3)(1/2) =", mp.nstr(W4M, 40))
# paper's Omega(1/2) Li4 combination:
C = (-21*mp.polylog(4,-2) - mp.mpf(65)/8*mp.polylog(4,mp.mpf(1)/3) + mp.mpf(13)/2*mp.polylog(4,mp.mpf(2)/3)
     + mp.mpf(137)/6*mp.polylog(4,mp.mpf(1)/2) + mp.mpf(39)/64*mp.polylog(4,mp.mpf(1)/9))
# predicted W4M = C - zeta3*Phi_wt1(1/2) + zeta2*Psi_wt2(1/2) = C + (3/8) zeta3 log3 + zeta2 Psi_wt2(1/2)
pred = C + mp.mpf(3)/8*mp.zeta(3)*mp.log(3) + mp.pi**2/6*Psi_wt2(c)
print("paper prediction  =", mp.nstr(pred, 40))
print("difference        =", mp.nstr(W4M - pred, 10))

# LLL: [W4M, Li4s, z3 l3, z3 l2, z2 Psi_wt2, products(26)]
l2, l3 = mp.log(2), mp.log(3)
z3 = mp.zeta(3); p2 = mp.pi**2
L2b = mp.polylog(2, mp.mpf(1)/3)
L3half = mp.polylog(3, mp.mpf(1)/2); L3b = mp.polylog(3, mp.mpf(1)/3); L3c = mp.polylog(3, mp.mpf(2)/3)
L2n = mp.polylog(2, mp.mpf(1)/9)
names = ["Li4(-2)", "Li4(1/3)", "Li4(2/3)", "Li4(1/2)", "Li4(1/9)", "z3 l3", "z3 l2", "z2 Psi2"]
xs = [mp.polylog(4,-2), mp.polylog(4,mp.mpf(1)/3), mp.polylog(4,mp.mpf(2)/3), mp.polylog(4,mp.mpf(1)/2),
      mp.polylog(4,mp.mpf(1)/9), z3*l3, z3*l2, p2/6*Psi_wt2(c)]
prods = [l2**4, l2**3*l3, l2**2*l3**2, l2*l3**3, l3**4,
         p2*l2**2, p2*l2*l3, p2*l3**2, mp.pi**4,
         z3*l2, z3*l3,
         L2b*l2**2, L2b*l2*l3, L2b*l3**2, L2b**2, p2*L2b,
         L3half*l2, L3b*l2, L3half*l3, L3b*l3, L3c*l2, L3c*l3,
         L2n*l2**2, L2n*l3**2, L2n**2]
pnames = (["l2^4","l2^3l3","l2^2l3^2","l2l3^3","l3^4","p2l2^2","p2l2l3","p2l3^2","pi4","z3l2p","z3l3p",
           "L2b l2^2","L2b l2l3","L2b l3^2","L2b^2","p2 L2b","L3half l2","L3b l2","L3half l3","L3b l3",
           "L3c l2","L3c l3","L2n l2^2","L2n l3^2","L2n^2"])
rels = find_relations([W4M] + xs + prods, 170, 8)
print("\nLLL shortest relations (a0 = coeff of W4M):")
for aa, resid in rels:
    mc = max(abs(x) for x in aa)
    print("   a0 =", aa[0], " maxcoef:", mc, " resid:", mp.nstr(resid, 3))
    if aa[0] != 0 and mc < 100000:
        coeff_Li4 = aa[1:6]
        print("     Li4 coeffs (x(-a0) sign conv):", coeff_Li4, " z3l3:", aa[6], " z3l2:", aa[7], " z2Psi2:", aa[8])
        print("     products:", [(pnames[i], aa[9+i]) for i in range(len(prods)) if aa[9+i]])
'''
open("<project>/research/08_theory/COR6X_independent_corner_symbols/step_06e.py","w").write(t)
import subprocess
r = subprocess.run(["nice","-n","10","<project>/.venv/bin/python","step_06e.py"],
                   capture_output=True, text=True, cwd="<project>/research/08_theory/COR6X_independent_corner_symbols", timeout=900)
print(r.stdout); print("ERR:", (r.stderr or "")[-1000:])
