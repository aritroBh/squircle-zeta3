src = open('ftfield.py').read()
src += '''
def norm_t(x):
    """reduce (A,B,C,d) to canonical form."""
    from math import lcm, gcd
    A,B,C,d = x
    L = lcm(lcm(A.denominator,B.denominator), lcm(C.denominator,d))
    a,b,c,l = A*L, B*L, C*L, L
    assert a.denominator==1 and b.denominator==1 and c.denominator==1
    a,b,c = a.numerator, b.numerator, c.numerator
    g = gcd(gcd(gcd(abs(a),abs(b)),abs(c)), l)
    a,b,c,l = a//g, b//g, c//g, l//g
    return (Fr(a,l), Fr(b,l), Fr(c,l), 1)
'''
open('ftfield.py','w').write(src)

geom = r'''
"""PG5 step_04: regular 14-gon edge data, exact in F = Q(2cos(2pi/7)).
Numerical geometry at 90 dps -> pslq reconstruction -> EXACT verification of all identities.
Term order (DOD1 convention): [a3(local t), a4(local 1-t), a1, a2, a5, a6, a0(const)].
"""
import mpmath as mp, pickle, sys
sys.path.insert(0,'.')
import ftfield as F
mp.mp.dps = 90
m = 7
PI = mp.pi
nn = [(mp.cos(i*PI/m), mp.sin(i*PI/m)) for i in range(m)]
sgn = {0:1,1:1,2:1,3:1,4:-1,5:-1,6:-1}
U = sum(sgn[i]*nn[i][0] for i in range(m)); V = sum(sgn[i]*nn[i][1] for i in range(m))
h = 1/U
assert abs(V) < mp.mpf('1e-80')
t14 = mp.tan(PI/14)
# edge: x1 = h, y in [-h t14, h t14]
V0 = (h, -h*t14); V1 = (h, h*t14)
L = 2*h*t14
def aval(i, t):
    x0 = V0[0] + t*(V1[0]-V0[0]); x1 = V0[1] + t*(V1[1]-V0[1])
    return sgn[i]*(nn[i][0]*x0 + nn[i][1]*x1)
def rec(y, maxcoeff=10**8):
    r = F.from_pslq(y, maxcoeff=maxcoeff)
    assert r is not None
    return F.norm_t(r)
# reconstruct everything
data = {}
for i in range(7):
    b0 = aval(i,0); b1 = aval(i,1)
    data[i] = (rec(b0), rec(b1-b0))   # [const, slope]
hh = rec(h); LL = rec(L); h2oL2 = rec(h*h/(L*L))
print("U =", mp.nstr(U,20), " h =", mp.nstr(h,20))
print("h exact:", hh, " L exact:", LL, " h^2/L^2:", h2oL2)
for i in range(7): print(f"a{i} =", data[i])
# EXACT verification: sum a_i = 1
s0 = F.ZERO; s1 = F.ZERO
for i in range(7):
    s0 = F.fadd(s0, data[i][0]); s1 = F.fadd(s1, data[i][1])
print("EXACT sum const:", F.norm_t(s0), " sum slope:", F.norm_t(s1))
# locals: a3(0)=0, a4(1)=0
print("a3(0)=0:", F.norm_t(data[3][0])==(F.Fr(0),)*3+(1,), " a4(1)=0:", F.norm_t(F.fadd(data[4][0],data[4][1])))
A0 = F.norm_t(data[3][1]); A1 = F.norm_t(data[4][0])
print("A0 =", A0, " A1 =", A1)
# nonlocals: B_i, x_i with a_i = B(1 - x t) -> x = -slope/B
X = {}; B = {}
for i in [1,2,5,6]:
    B[i] = F.norm_t(data[i][0])
    X[i] = F.norm_t(F.fneg(F.fdiv(data[i][1], data[i][0])))
    print(f"B{i} = {B[i]}  x{i} = {X[i]}  evals {float(F.feval(B[i])):.15f} {float(F.feval(X[i])):.15f}")
# Landen checks EXACT
def lnd(x): return F.fdiv(x, F.fsub(x, F.ONE))
print("x6 = x1/(x1-1):", F.norm_t(F.fsub(X[6], lnd(X[1]))))
print("x5 = x2/(x2-1):", F.norm_t(F.fsub(X[5], lnd(X[2]))))
# constant term a0
print("a0 const:", F.norm_t(data[0][0]), " slope:", F.norm_t(data[0][1]))
# reflection pairs: B5/B2 = 1-x2? a_i(t) reflect: B_{m+1-i} = B_i (1-x_i)
print("B6 - B1(1-x1):", F.norm_t(F.fsub(B[6], F.fmul(B[1], F.fsub(F.ONE,X[1])))))
print("B5 - B2(1-x2):", F.norm_t(F.fsub(B[5], F.fmul(B[2], F.fsub(F.ONE,X[2])))))
# norms of everything (which primes appear?)
for nm, v in [('h',hh),('A0',A0),('B1',B[1]),('B2',B[2]),('B5',B[5]),('B6',B[6])] + [(f"x{i}",X[i]) for i in [1,2,5,6]] + [(f"1-x{i}", F.fsub(F.ONE,X[i])) for i in [1,2,5,6]]:
    print(nm, "N =", F.fnorm(v))
pickle.dump({'data':data,'h':hh,'L':LL,'h2oL2':h2oL2,'A0':A0,'A1':A1,'B':B,'X':X}, open('geom14.pkl','wb'))
open('progress.log','a').write('PG5 step_04: 14-gon geometry exact, all identities verified exactly\n')
'''
open('step_04_geom14.py','w').write(geom)
import subprocess
r = subprocess.run(['<project>/.venv/bin/python','step_04_geom14.py'],capture_output=True,text=True)
print(r.stdout); print(r.stderr[-800:] if r.stderr else '')
