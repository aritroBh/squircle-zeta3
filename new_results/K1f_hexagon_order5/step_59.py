code = r'''
# step_27.py: local Beta table Jloc(N;a,b) = int_0^1 t^N log^a(t) log^b(1-t) dt, a+b<=5.
# Symbolic in N via harmonic symbols Hs(r) = H_{N+1}^{(r)}, T2(r) = 1/(N+2)^r, den = 1/(N+1).
# Constants z2,z3,z4=2z2^2/5,z5. Verified vs mpmath at 50 digits for random N.
import sympy as sp, pickle
WD = '<project>/research/08_theory/K1f_hexagon_order5'
u, v = sp.symbols('u v')
z2, z3, z5 = sp.symbols('z2 z3 z5')
z4 = sp.Rational(2,5)*z2**2
Hs = sp.symbols('H1 H2 H3 H4 H5')          # H_{N+1}^{(r)}, r=1..5
T2 = sp.symbols('T1 T2 T3 T4 T5')          # 1/(N+2)^r
DEN = sp.Symbol('DEN')                     # 1/(N+1)

def lg_series_at_N1():
    # lnGamma(N+1+z)-lnGamma(N+1) = sum_j c_j z^j ; H_{N+1}^{(r)} = Hs[r-1]
    return [Hs[0], (z2 - Hs[1])/2, (-2*z3 + 2*Hs[2])/6, (6*z4 - 6*Hs[3])/24, (-24*z5 + 24*Hs[4])/120]
def lg_series_at_1():
    return [sp.Integer(0), z2/2, -z3/3, z4/4, -z5/5]
def lg_series_at_N2():
    # H_{N+2}^{(r)} = Hs[r-1] + T2[r-1]
    H = [Hs[r] + T2[r] for r in range(5)]
    return [H[0], (z2 - H[1])/2, (-2*z3 + 2*H[2])/6, (6*z4 - 6*H[3])/24, (-24*z5 + 24*H[4])/120]

A = lg_series_at_N1(); B = lg_series_at_1(); C = lg_series_at_N2()
S = sp.expand(sum(A[j]*u**(j+1) + B[j]*v**(j+1) - C[j]*(u+v)**(j+1) for j in range(5)))
E = sp.Integer(1); P = sp.Integer(1)
for j in range(1, 6):
    P = sp.expand(P*S)
    E = sp.expand(E + P/sp.factorial(j))
tab = {}
for a in range(6):
    for b in range(6-a):
        co = sp.expand(E.coeff(u, a).coeff(v, b))
        tab[(a, b)] = sp.expand(DEN*sp.factorial(a)*sp.factorial(b)*co)
# gamma should not appear (it doesn't: lg_series has no gamma since we removed it: check H1 coefficient cancels)
print("built Jloc table, entries:", len(tab))
for key in [(1,0),(0,1),(2,0),(1,1),(3,2),(4,1),(5,0),(0,5)]:
    print(f"Jloc(N;{key}) =", sp.collect(tab[key], [z5,z2*z3,z3,z2]))
pickle.dump({k: str(v) for k,v in tab.items()}, open(f'{WD}/jloc_sym.pkl','wb'))
print("saved jloc_sym.pkl")

# numeric verification
import mpmath as mp
mp.mp.dps = 50
zv = {z2: sp.pi**2/6, z3: sp.zeta(3), z5: sp.zeta(5)}
import random
random.seed(1)
worst = mp.mpf(0)
for trial in range(12):
    N = random.randint(0, 40)
    subs = {Hs[r-1]: sp.harmonic(N+1, r) for r in range(1,6)}
    subs.update({T2[r-1]: sp.Rational(1,(N+2)**r) for r in range(1,6)})
    subs[DEN] = sp.Rational(1, N+1)
    subs.update(zv)
    for (a,b), e in tab.items():
        val = mp.mpf(sp.N(e.subs(subs, simultaneous=True), 50))
        num = mp.quad(lambda t: t**N*mp.log(t)**a*mp.log(1-t)**b, [0, mp.mpf('0.5'), 1])
        d = abs(val-num); worst = max(worst, d)
print("verification: 12 random N x 21 (a,b), worst |diff| =", mp.nstr(worst, 3))
'''
wd = '<project>/research/08_theory/K1f_hexagon_order5'
open(f'{wd}/step_27.py','w').write(code)
import subprocess
PY='<project>/.venv/bin/python'
r = subprocess.run([PY, f'{wd}/step_27.py'], capture_output=True, text=True, cwd=wd, timeout=900)
print(r.stdout)
print(r.stderr[-1500:] if r.returncode else '')
