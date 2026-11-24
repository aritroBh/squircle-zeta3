import os
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
src = '''"""oct5d3_atoms.py — numerical evaluation of all atoms at working precision.

Letters: lam = log 2, L = log(1+sqrt2). z_k = zeta(k), Lchi_k = L(chi8,k).
Li32(x,1) = sum_{n>m>0} x^n/(n^3 m^2) = sum_n x^n H^{(2)}_{n-1}/n^3 (single series).
"""
import mpmath as mp

def lchi(k):
    return sum(s*mp.zeta(k, mp.mpf(a)/8) for a, s in ((1, 1), (7, 1), (3, -1), (5, -1)))/mp.mpf(8)**k

def li32_x1(x):
    """Li_{3,2}(x,1) = sum_{n>=2} x^n H^{(2)}_{n-1}/n^3, |x|<1."""
    tot = mp.mpf(0); H2 = mp.mpf(0); xp = x; n = 1
    thresh = mp.mpf(10)**(-mp.mp.dps-15)
    while True:
        H2 += mp.mpf(1)/(n*n)
        n += 1
        xp *= x
        tot += xp*H2/mp.mpf(n)**3
        if n > 10 and abs(xp) < thresh: break
    return tot

def base_atoms():
    r2 = mp.sqrt(2)
    A = {}
    A['r2'] = r2
    A['lam'] = mp.log(2); A['L'] = mp.log(1+r2)
    A['pi2'] = mp.pi**2
    for k in (2, 3, 4, 5):
        A[f'z{k}'] = mp.zeta(k)
    A['Lchi3'] = lchi(3); A['Lchi5'] = lchi(5)
    b = r2-1
    A['Li2(b)'] = mp.polylog(2, b)
    for nm, z in [('b', b), ('-b', -b), ('1/sqrt2', 1/r2)]:
        A[f'Li3({nm})'] = mp.polylog(3, z)
    for nm, z in [('b', b), ('-b', -b), ('b^2', b*b), ('1/sqrt2', 1/r2), ('1/2', mp.mpf('0.5')),
                  ('1-1/sqrt2', 1-1/r2), ('2-sqrt2', 2-r2)]:
        A[f'Li4({nm})'] = mp.polylog(4, z)
    for nm, z in [('b', b), ('-b', -b), ('-b^2', -b*b), ('1/sqrt2', 1/r2), ('1/2', mp.mpf('0.5')),
                  ('1-1/sqrt2', 1-1/r2), ('2-sqrt2', 2-r2),
                  ('2b', 2*b), ('(2-sqrt2)/4', (2-r2)/4), ('(2+sqrt2)/4', (2+r2)/4),
                  ('-b/2', -b/2), ('(4-3sqrt2)/8', mp.mpf('0.5')-3*r2/8),
                  ('b^4', b**4), ('1-b^4', 1-b**4)]:
        A[f'Li5({nm})'] = mp.polylog(5, z)
    A['Li32(b,1)'] = li32_x1(b)
    A['Li32(1/sqrt2,1)'] = li32_x1(1/r2)
    return A

def prod(A, spec):
    """spec: list of (name, power) -> product."""
    v = mp.mpf(1)
    for nm, pw in spec:
        v *= A[nm]**pw
    return v

# planner's 22 invisible weight-5 atoms: (label, spec)
INVISIBLE5 = [
    ("z2*Li3(b)", [('z2',1),('Li3(b)',1)]),
    ("z2*Li3(-b)", [('z2',1),('Li3(-b)',1)]),
    ("z2*Li3(1/sqrt2)", [('z2',1),('Li3(1/sqrt2)',1)]),
    ("z2*z3", [('z2',1),('z3',1)]),
    ("z2*Lchi3", [('z2',1),('Lchi3',1)]),
    ("z2*lam^3", [('z2',1),('lam',3)]),
    ("z2*lam^2*L", [('z2',1),('lam',2),('L',1)]),
    ("z2*lam*L^2", [('z2',1),('lam',1),('L',2)]),
    ("z2*L^3", [('z2',1),('L',3)]),
    ("z2*Li2(b)*lam", [('z2',1),('Li2(b)',1),('lam',1)]),
    ("z2*Li2(b)*L", [('z2',1),('Li2(b)',1),('L',1)]),
    ("z2*pi2*lam", [('z2',1),('pi2',1),('lam',1)]),
    ("z2*pi2*L", [('z2',1),('pi2',1),('L',1)]),
    ("z3*pi2", [('z3',1),('pi2',1)]),
    ("z3*lam^2", [('z3',1),('lam',2)]),
    ("z3*lam*L", [('z3',1),('lam',1),('L',1)]),
    ("z3*L^2", [('z3',1),('L',2)]),
    ("z3*Li2(b)", [('z3',1),('Li2(b)',1)]),
    ("z4*lam", [('z4',1),('lam',1)]),
    ("z4*L", [('z4',1),('L',1)]),
    ("z5", [('z5',1)]),
    ("Lchi5", [('Lchi5',1)]),
]

# lower-weight bases
W4_BASIS = (
    [(f"Li4({nm})", [(f'Li4({nm})',1)]) for nm in ['b','-b','b^2','1/sqrt2','1/2','1-1/sqrt2','2-sqrt2']]
    + [(f"lam^{i}*L^{4-i}", [('lam',i),('L',4-i)]) for i in range(5)]
    + [(f"pi2*lam^{i}*L^{2-i}", [('pi2',1),('lam',i),('L',2-i)]) for i in range(3)]
    + [("pi2^2", [('pi2',2)]), ("Li2(b)^2", [('Li2(b)',2)])]
    + [(f"Li2(b)*lam^{i}*L^{2-i}", [('Li2(b)',1),('lam',i),('L',2-i)]) for i in range(3)]
    + [("pi2*Li2(b)", [('pi2',1),('Li2(b)',1)])]
    + [("z3*lam", [('z3',1),('lam',1)]), ("z3*L", [('z3',1),('L',1)]),
       ("Lchi3*lam", [('Lchi3',1),('lam',1)]), ("Lchi3*L", [('Lchi3',1),('L',1)])]
    + [(f"Li3({nm})*{lg}", [(f'Li3({nm})',1),(lg,1)]) for nm in ['b','-b','1/sqrt2'] for lg in ['lam','L']]
)
W3_BASIS = (
    [(f"Li3({nm})", [(f'Li3({nm})',1)]) for nm in ['b','-b','1/sqrt2']]
    + [("z3", [('z3',1)]), ("Lchi3", [('Lchi3',1)])]
    + [(f"lam^{i}*L^{3-i}", [('lam',i),('L',3-i)]) for i in range(4)]
    + [(f"Li2(b)*{lg}", [('Li2(b)',1),(lg,1)]) for lg in ['lam','L']]
    + [(f"pi2*{lg}", [('pi2',1),(lg,1)]) for lg in ['lam','L']]
)
W2_BASIS = [("pi2", [('pi2',1)]), ("lam^2", [('lam',2)]), ("lam*L", [('lam',1),('L',1)]),
            ("L^2", [('L',2)]), ("Li2(b)", [('Li2(b)',1)])]
W1_BASIS = [("lam", [('lam',1)]), ("L", [('L',1)])]

def eval_table(dps):
    mp.mp.dps = dps
    return base_atoms()

if __name__ == "__main__":
    import sys, time
    dps = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    mp.mp.dps = dps
    t0 = time.time()
    A = base_atoms()
    print("base atoms evaluated in %.1fs" % (time.time()-t0))
    print("Li32(b,1)      =", mp.nstr(A['Li32(b,1)'], 40))
    print("Li32(1/sqrt2,1)=", mp.nstr(A['Li32(1/sqrt2,1)'], 40))
    print("Lchi5 =", mp.nstr(A['Lchi5'], 40))
    print("Lchi3 =", mp.nstr(A['Lchi3'], 40))
    # cross-check Li32 against OCT5D2 430-digit values
    import pickle
    av = pickle.load(open("<project>/research/08_theory/OCT5D2_octagon_depth_test/atom_vals430.pkl","rb"))
    print("Li32(b,1) - OCT5D2  =", mp.nstr(abs(A['Li32(b,1)']-mp.mpf(av['L32_b1'])), 3))
    print("Li32(1/sqrt2,1)-OCT5D2=", mp.nstr(abs(A['Li32(1/sqrt2,1)']-mp.mpf(av['L32_h1'])), 3))
    # Lchi5 vs M4e 1650-digit value
    import json
    l5 = json.load(open("<project>/research/08_theory/M4f_w5_extras/Lchi8_5.json"))
    print("Lchi8_5.json keys:", list(l5.keys())[:5])
'''
open("oct5d3_atoms.py","w").write(src)
import subprocess, sys
r = subprocess.run([sys.executable, "oct5d3_atoms.py", "120"], capture_output=True, text=True)
print(r.stdout); print(r.stderr[-1500:] if r.returncode else "OK")
