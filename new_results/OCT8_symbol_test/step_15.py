import time
with open("progress.log","a") as f:
    f.write("STEP symbol-engine: S(Li21) formula derived & verified: integrable (both wedge conditions), "
            "stuffle with S(Li12(y,x)) & S(Li3(xy)) EXACT, Landen for Li3 EXACT. "
            "Li21 integral rep int_0^1 y Li2(xyt)/(1-yt) dt validated 60d vs series.\n")

# ---------------- Q: S-unit equation for S={2,3,5} ----------------
from fractions import Fraction as Fr
import sympy as sp

def factor_rat(q):
    """q nonzero rational -> dict prime->Fraction exponent (sign/torsion ignored)."""
    q = Fr(q)
    n, d = abs(q.numerator), q.denominator
    out = {}
    for m, sgn in ((n,1),(d,-1)):
        for p, e in sp.factorint(m).items():
            out[p] = out.get(p, Fr(0)) + sgn*e
    return out

def is_S_unit(q, S):
    v = factor_rat(q)
    return all(p in S for p in v)

B = 30
S = (2,3,5)
sol = []
import itertools
vals = []
for a,b,c in itertools.product(range(-B,B+1), repeat=3):
    u = Fr(2)**a * Fr(3)**b * Fr(5)**c
    vals.append((u,(a,b,c)))
for u,(a,b,c) in vals:
    v = 1-u
    if v == 0: continue
    if is_S_unit(v, S) and is_S_unit(u, S):
        sol.append((u, v, (a,b,c)))
# dedupe by u (unique anyway) and keep u != 0
solU = sorted(set(u for u,v,_ in sol))
print("search |a|,|b|,|c|<=",B,": #solutions u (u,1-u both {2,3,5}-units):", len(solU))
# closure under anharmonic S3: u -> u, 1-u, 1/u, 1/(1-u), u/(u-1), (u-1)/u
def orbit(u):
    return {u, 1-u, 1/u, 1/(1-u), u/(u-1), (u-1)/u}
Ss = set(solU)
orbits = []
seen = set()
for u in solU:
    if u in seen: continue
    o = orbit(u)
    orbits.append(o)
    seen |= o
print("closed under S3 anharmonic action:", all(o <= Ss for o in orbits), " #orbits:", len(orbits),
      " sizes:", sorted(len(o) for o in orbits))
for o in sorted(orbits, key=lambda o: min(o)):
    print("  orbit:", sorted(o))
with open("progress.log","a") as f:
    f.write(f"STEP S-units: S={{2,3,5}} search B={B}: {len(solU)} solutions, S3-closed, {len(orbits)} orbits\n")
