code = r'''
# step_29.py: build the weight-5 non-local J-table:
#   J(n;a,b,c) = int_0^1 t^n X^a Y^b W^c dt = c!(-1)^c sum_{m>=c} x^m e_{c-1}(m-1)/m * Jloc(n+m;a,b)
# with x = -kappa symbolic, atoms in the verified HPL basis.
import sympy as sp, pickle, sys, time
WD = '<project>/research/08_theory/K1f_hexagon_order5'
sys.path.insert(0, WD)
import esum5
esum5.load_atom_table(WD)
x = esum5.x; Hf = esum5.Hf
m = sp.Symbol('m')
z2, z3, z5 = sp.symbols('z2 z3 z5')

# Jloc symbolic table
jloc_raw = pickle.load(open(f'{WD}/jloc_sym.pkl','rb'))
Hs = sp.symbols('H1 H2 H3 H4 H5')
T2 = sp.symbols('T1 T2 T3 T4 T5')
DEN = sp.Symbol('DEN')
jloc = {k: sp.sympify(v) for k, v in jloc_raw.items()}

def jloc_sub(n, a, b):
    """Jloc(n+m; a, b) as expr in summation index m."""
    e = jloc[(a, b)]
    subs = {DEN: 1/(m+n+1)}
    for r in range(1, 6):
        subs[Hs[r-1]] = Hf(r, m+n+1)
        subs[T2[r-1]] = 1/(m+n+2)**r
    return sp.expand(e.subs(subs, simultaneous=True))

def esym(c, mv):
    """e_{c-1}(1,1/2,...,1/(mv-1)) with Hf(j, mv-1)."""
    H = [None] + [Hf(j, mv-1) for j in range(1, 5)]
    if c == 1: return sp.Integer(1)
    if c == 2: return H[1]
    if c == 3: return (H[1]**2 - H[2])/2
    if c == 4: return (H[1]**3 - 3*H[1]*H[2] + 2*H[3])/6
    if c == 5: return (H[1]**4 - 6*H[1]**2*H[2] + 3*H[2]**2 + 8*H[1]*H[3] - 6*H[4])/24
    raise ValueError(c)

def J(n, a, b, c):
    g = sp.expand(esym(c, m)/m * jloc_sub(n, a, b))
    return sp.expand(sp.factorial(c)*(-1)**c * esum5.esum(g, c, m))

if __name__ == '__main__':
    t0 = time.time()
    for (n,a,b,c) in [(0,0,0,1),(2,1,1,2),(0,0,0,2)]:
        Je = J(n,a,b,c)
        print(f"J({n};{a},{b},{c}) built [{time.time()-t0:.0f}s], terms={len(sp.Add.make_args(Je))}")
        print("   ", Je)
    pickle.dump({'test': 1}, open(f'{WD}/jtab_test.pkl','wb'))
'''
wd = '<project>/research/08_theory/K1f_hexagon_order5'
open(f'{wd}/step_29.py','w').write(code)
import subprocess
PY='<project>/.venv/bin/python'
r = subprocess.run([PY, f'{wd}/step_29.py'], capture_output=True, text=True, cwd=wd, timeout=1700)
print(r.stdout[-4000:])
print(r.stderr[-2500:] if r.returncode else '')
