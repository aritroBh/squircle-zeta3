
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, ".")
import mpmath as mp
mp.mp.dps = 60
import reflib as R
# weight-2 part of Psi(1/2) from the G1 fit (1,0,0): terms with |w|=2
import json, re
G1 = "<project>/research/10_frontier/G1_zeta3_any_angle"
D = json.load(open(G1 + "/hpl_fit_k3_(1,0,0).json"))
psi2 = mp.mpf(0)
for key, val in D:
    m = re.match(r"c\^(-?\d+) H\[(.*)\]", key)
    if not m: continue
    e = int(m.group(1))
    w = tuple(int(z) for z in m.group(2).split(",")) if m.group(2).strip() else ()
    if len(w) != 2: continue
    psi2 += mp.mpf(str(val)) * mp.mpf(1)/2**e * R.hpl_value(w, mp.mpf(1)/2) / (1 if e >= 0 else 1)
    print(key, val, "->", mp.nstr(R.hpl_value(w, mp.mpf(1)/2), 8))
print("Psi2(1/2) mine:", mp.nstr(psi2, 25))
ref = mp.mpf(open("../COR5_corner_classicality/hex_Psi2.txt").read().strip())
print("hex_Psi2.txt    :", mp.nstr(ref, 25))
print("|diff| =", mp.nstr(abs(psi2-ref), 3))
# full Psi(1/2) check too
psi = mp.mpf(0)
for key, val in D:
    m = re.match(r"c\^(-?\d+) H\[(.*)\]", key)
    if key == "1": e, w = 0, ()
    elif key == "c": e, w = 1, ()
    else:
        e = int(m.group(1)); w = tuple(int(z) for z in m.group(2).split(",")) if m.group(2).strip() else ()
    v = mp.mpf(1) if not w else R.hpl_value(w, mp.mpf(1)/2)
    psi += mp.mpf(str(val)) * mp.mpf(1)/2**e * v
print("Psi(1/2) mine:", mp.nstr(psi, 25), "vs hex_Psi1:", open("../COR5_corner_classicality/hex_Psi1.txt").read().strip()[:27])
print("also closed form -(5/2)log2 + (3/4)log3:", mp.nstr(-mp.mpf(5)/2*mp.log(2)+mp.mpf(3)/4*mp.log(3), 25))
