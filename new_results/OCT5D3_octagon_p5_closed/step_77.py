import os
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
src = open("oct5d3_atoms.py").read()
# replace INVISIBLE5 with the 21-atom version (12 weight-3, no Lchi3) + a 22-atom variant
i0 = src.find("# planner's 22 invisible weight-5 atoms")
i1 = src.find("# lower-weight bases")
new = '''# planner's invisible weight-5 basis (21 atoms: zeta2 x 12 weight-3 [paper convention,
# chi3(b) spans the L(chi8,3) direction], zeta3 x 5, zeta4 x 2, zeta5, Lchi5)
W3_12 = [('Li3(b)',1), ('Li3(-b)',1), ('Li3(1/sqrt2)',1), ('z3',1),
         ('lam',3), ('lam',2), ('lam',1), ('L',3), ('Li2(b)',1), ('Li2(b)',1), ('pi2',1), ('pi2',1)]
INVISIBLE5 = [
    ("z2*Li3(b)", [('z2',1),('Li3(b)',1)]),
    ("z2*Li3(-b)", [('z2',1),('Li3(-b)',1)]),
    ("z2*Li3(1/sqrt2)", [('z2',1),('Li3(1/sqrt2)',1)]),
    ("z2*z3", [('z2',1),('z3',1)]),
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
INVISIBLE5_PLUS = INVISIBLE5 + [("z2*Lchi3", [('z2',1),('Lchi3',1)])]

'''
src = src[:i0] + new + src[i1:]
open("oct5d3_atoms.py","w").write(src)
print("updated")
import subprocess, sys
r = subprocess.run([sys.executable, "-c", "import sys; sys.path.insert(0,'.'); from oct5d3_atoms import INVISIBLE5, INVISIBLE5_PLUS; print(len(INVISIBLE5), len(INVISIBLE5_PLUS))"], capture_output=True, text=True)
print(r.stdout, r.stderr[-300:] if r.returncode else "")
