import os, json, shutil, re
from fractions import Fraction as Fr
CWD = "<project>/research/08_theory/G2R_allN_proof_referee"
RB = os.path.join(CWD, "run_base")
FILES = [f for f in os.listdir(RB) if f != "mu4_zsum.pkl"]

def fresh(dst):
    os.makedirs(dst, exist_ok=True)
    for fn in FILES:
        shutil.copy2(os.path.join(RB, fn), os.path.join(dst, fn))

log = []

# ---- (a) one rational coefficient of one weight-5 term of R4 ----
dst = os.path.join(CWD, "mut_a"); fresh(dst)
p = os.path.join(dst, "hpl_fit_k4_(0,0,0).json")
d = json.load(open(p))
done = False
for i, (nm, v) in enumerate(d):
    m = re.match(r'c\^(-?\d+) He\[(.*)\]', nm)
    if not m: continue
    w = tuple(int(t) for t in m.group(2).split(','))
    if len(w) == 5:
        newv = str(Fr(v) + 1)
        log.append(f"(a) R4 weight-5 term {nm}: coeff {v} -> {newv}")
        d[i] = [nm, newv]; done = True
        break
assert done
json.dump(d, open(p, "w"))

# ---- (b) one coefficient of Psi2 ----
dst = os.path.join(CWD, "mut_b"); fresh(dst)
p = os.path.join(dst, "hpl_fit_k4_(1,0,0).json")
d = json.load(open(p))
done = False
for i, (nm, v) in enumerate(d):
    if nm == 'c^0 He[1]':
        newv = str(Fr(v) + Fr(1, 3))
        log.append(f"(b) Psi2 term {nm}: coeff {v} -> {newv}")
        d[i] = [nm, newv]; done = True
        break
assert done
json.dump(d, open(p, "w"))

# ---- (c) one coefficient of F1 ----
dst = os.path.join(CWD, "mut_c"); fresh(dst)
p = os.path.join(dst, "b4_pi4_part.json")
d = json.load(open(p))
done = False
for i, (cof, kind) in enumerate(d["F1"]):
    if cof == "7*c/16":
        log.append(f"(c) F1 artanh coeff {cof} -> 5*c/16")
        d["F1"][i] = ["5*c/16", kind]; done = True
        break
assert done
json.dump(d, open(p, "w"))

# ---- (d) the prefactor exponent e of one term ----
dst = os.path.join(CWD, "mut_d"); fresh(dst)
p = os.path.join(dst, "hpl_fit_k4_(0,0,0).json")
d = json.load(open(p))
done = False
for i, (nm, v) in enumerate(d):
    if nm == 'c^-2 He[1, 0, 1]':
        newnm = 'c^-3 He[1, 0, 1]'
        log.append(f"(d) R4 term exponent: {nm} -> {newnm} (coeff {v} unchanged)")
        d[i] = [newnm, v]; done = True
        break
assert done
json.dump(d, open(p, "w"))

for line in log: print(line)
with open(CWD + "/progress.log", "a") as fh:
    fh.write("step_mutants_setup: " + " | ".join(log) + "\n")
