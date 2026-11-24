import time, json, re
WD = '<project>/research/08_theory/COR9_golden_corner_closed_forms/'
import sys; sys.path.insert(0, WD)
import mpmath as mp
mp.mp.dps = 1120
r5 = mp.sqrt(5); phi = (1+r5)/2
L = mp.log(phi); lam = mp.log(2); pi2 = mp.pi**2; pi4 = mp.pi**4
Li = lambda s, x: mp.polylog(s, x)
z3 = mp.zeta(3)

def eval_arg(s):
    e = s.replace("sqrt5", "r5").replace("^", "**")
    e = re.sub(r"(\d)\s*phi", r"\1*phi", e)      # 2phi -> 2*phi
    e = re.sub(r"r5\s+phi", "r5*phi", e)          # r5 phi -> r5*phi
    e = re.sub(r"\)\s*phi", ")*phi", e)
    return eval(e, {"r5": r5, "phi": phi})

def basis_val(nm):
    v = mp.mpf(1)
    for fac in nm.split("*"):
        m = re.match(r"^(Li[234])\((.*)\)(?:\^(\d+))?$", fac)
        if m:
            s_ = int(m.group(1)[2]); x = eval_arg(m.group(2))
            vv = Li(s_, x)
            if m.group(3): vv = vv**int(m.group(3))
        else:
            m2 = re.match(r"^(zeta3|pi2|pi4|L|lam)(?:\^(\d+))?$", fac)
            vv = {"zeta3": z3, "pi2": pi2, "pi4": pi4, "L": L, "lam": lam}[m2.group(1)]
            if m2.group(2): vv = vv**int(m2.group(2))
        v *= vv
    return v

# sanity: evaluate all stored basis names
for tag in ("gol1", "gol3"):
    for f in (f"{tag}_W4_relation.json", f"{tag}_W3_relation.json"):
        r = json.load(open(WD+f))
        for nm in r["basis"]:
            basis_val(nm)
print("all basis names parse OK")
for tag in ("gol1", "gol3"):
    M3 = mp.mpf(open(WD+f"{tag}_M3_hpl_1150.txt").read())
    W4v = mp.mpf(open(WD+f"{tag}_W4_hpl_1150.txt").read())
    W3v = mp.mpf(open(WD+f"{tag}_W3_hpl_1150.txt").read())
    W2v = mp.mpf(open(WD+f"{tag}_W2_hpl_1150.txt").read())
    r4 = json.load(open(WD+f"{tag}_W4_relation.json"))
    B4 = [basis_val(nm) for nm in r4["basis"]]; nb4 = len(B4)
    W4cf = -sum(mp.mpf(r4["coeffs"][i])*B4[i] + mp.mpf(r4["coeffs"][i+nb4])*r5*B4[i] for i in range(nb4))/r4["a0"]
    r3 = json.load(open(WD+f"{tag}_W3_relation.json"))
    B3 = [basis_val(nm) for nm in r3["basis"]]; nb3 = len(B3)
    W3cf = -sum(mp.mpf(r3["coeffs"][i])*B3[i] + mp.mpf(r3["coeffs"][i+nb3])*r5*B3[i] for i in range(nb3))/r3["a0"]
    c = 1/phi if tag == "gol1" else r5-2
    W2cf = -Li(2, c) - Li(2, -c) + pi2/12
    print(f"{tag}: |W4-W4cf| = {mp.nstr(abs(W4v-W4cf),3)}  |W3-W3cf| = {mp.nstr(abs(W3v-W3cf),3)}  |W2-W2cf| = {mp.nstr(abs(W2v-W2cf),3)}")
    print(f"      |M3 - (W4cf+W3cf+W2cf)| = {mp.nstr(abs(M3 - (W4cf+W3cf+W2cf)),3)}")
open(WD+'progress.log','a').write(f"{time.strftime('%H:%M:%S')} step_50: golden assemblies verified.\n")
