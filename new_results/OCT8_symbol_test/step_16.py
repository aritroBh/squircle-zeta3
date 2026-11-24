from fractions import Fraction as Fr
import itertools, time

def strip235(n):
    n = abs(n)
    for p in (2,3,5):
        while n % p == 0:
            n //= p
    return n

B = 14
sol = set()
t0 = time.time()
for a,b,c in itertools.product(range(-B,B+1), repeat=3):
    u = Fr(2)**a * Fr(3)**b * Fr(5)**c
    v = 1-u
    if v == 0: continue
    # v is S-unit iff numerator and denominator of v are {2,3,5}-smooth
    if strip235(v.numerator) == 1 and strip235(v.denominator) == 1:
        sol.add(u)
print("B =", B, " time %.1fs"%(time.time()-t0), " #solutions:", len(sol))

def orbit(u):
    return {u, 1-u, 1/u, 1/(1-u), u/(u-1), (u-1)/u}
Ss = set(sol)
seen = set(); orbits = []
for u in sorted(Ss):
    if u in seen: continue
    o = orbit(u); orbits.append(o); seen |= o
closed = all(o <= Ss for o in orbits)
print("S3-closed:", closed, " #orbits:", len(orbits))
for o in sorted(orbits, key=lambda o: min(o)):
    rep = sorted(o)
    print("  size %d:"%len(o), rep)
with open("progress.log","a") as f:
    f.write(f"STEP S-units redo: B={B}, {{2,3,5}}-unit eqn u+v=1: {len(sol)} sols, S3-closed={closed}, orbits={len(orbits)}\n")
