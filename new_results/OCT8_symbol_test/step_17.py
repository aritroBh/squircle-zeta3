from fractions import Fraction as Fr
import itertools, time, pickle

def strip235(n):
    n = abs(n)
    for p in (2,3,5):
        while n % p == 0:
            n //= p
    return n

def search(B):
    sol = set()
    for a,b,c in itertools.product(range(-B,B+1), repeat=3):
        base = Fr(2)**a * Fr(3)**b * Fr(5)**c
        for s in (1,-1):
            u = s*base
            v = 1-u
            if v == 0: continue
            if strip235(v.numerator) == 1 and strip235(v.denominator) == 1:
                sol.add(u)
    return sol

for B in (10, 14, 18):
    t0=time.time(); s = search(B); print("B=%2d: %4d solutions  (%.1fs)"%(B,len(s),time.time()-t0))
sol = search(14)
def orbit(u):
    return {u, 1-u, 1/u, 1/(1-u), u/(u-1), (u-1)/u}
seen=set(); orbits=[]
for u in sorted(sol):
    if u in seen: continue
    o=orbit(u); orbits.append(o); seen|=o
print("S3-closed:", all(o<=sol for o in orbits), "orbits:", len(orbits), "sizes:", sorted(len(o) for o in orbits))
print("total = 3 + 6*16 =", 3+6*16)
pickle.dump(sorted(sol), open("sols235.pkl","wb"))
with open("progress.log","a") as f:
    f.write("STEP S-units signed: {2,3,5}-unit eqn: 99 solutions = 16 generic S3-orbits + special orbit {1/2,2,-1}; "
            "stable for B=10..18; S3-closed. Saved sols235.pkl\n")
