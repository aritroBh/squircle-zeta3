import os
work = "<project>/research/08_theory/AV5b_layer_theorem"
os.chdir(work)
ver = r'''
"""verify_av5b.py - AV5b layer theorem verification suite (<10 min)."""
import mpmath as mp
import time

t_start = time.time()
mp.mp.dps = 50
z3, z2, z4 = mp.zeta(3), mp.zeta(2), mp.zeta(4)
pi2 = mp.pi**2
results = []
def check(name, val, target, tol, label):
    ok = abs(val - target) < mp.mpf(tol)
    results.append((name, ok, val, target, label))
    print("[%s] %s: |val-target|=%.2e (tol %s)  %s" % ("PASS" if ok else "FAIL", name, abs(val-target), tol, label), flush=True)

# ---------- Task A: softplus moment formulas ----------
g1 = lambda s: mp.log(1+mp.exp(-2*s))
for j in [0,1,2,3]:
    v = mp.quad(lambda s: s**j*g1(s), [0, mp.inf])
    t = mp.factorial(j)*mp.eta(j+2)/2**(j+1)
    check("A.moment j=%d" % j, v, t, '1e-40', "PROVED+NUMERICAL(40d)")
v = mp.quad(lambda s: g1(s)**2, [0, mp.inf])
check("A.int g1^2 = zeta(3)/8", v, z3/8, '1e-40', "PROVED+NUMERICAL(40d)")
S = mp.mpf(0)
N = 4000
for jj in range(1, N+1):
    S += mp.fsum((-1)**(jj+k)/(jj*k*(jj+k)**2) for k in range(1, N+1))
I_num = S/4
I_dir = mp.quad(lambda s: s*g1(s)**2, [0, mp.inf])
check("A.int s g1^2 = coloured-MZV series", I_num, I_dir, '1e-6', "PROVED(representation)")
J = mp.quad(lambda s: s*g1(s)**2/2, [0, mp.inf]) + mp.quad(lambda s: g1(s)**3/6, [0, mp.inf])
check("A.collapse J = zeta(4)/32", J, z4/32, '1e-40', "NUMERICAL(40d)")

# ---------- Task B: squircle (m=2) ----------
def P_squircle(p):
    q = p/(p-1)
    return 8*mp.quad(lambda th: (mp.cos(th)**q + mp.sin(th)**q)**(1/q), [0, mp.pi/4])
def neville(pts):
    n = len(pts); Q = [[mp.mpf(0)]*n for _ in range(n)]
    for i in range(n): Q[i][0] = pts[i][1]
    for j in range(1,n):
        for i in range(n-j):
            xi, xj = pts[i][0], pts[i+j][0]
            Q[i][j] = (xj*Q[i][j-1] - xi*Q[i+1][j-1])/(xj-xi)
    return Q[0][n-1]
taus = [mp.mpf('0.03')/2**k for k in range(6)]
Ps = [P_squircle(1/(1-t)) - 4*mp.sqrt(2) for t in taus]
c2 = neville([(t, v/t**2) for t,v in zip(taus,Ps)])
rem = [(t, v - c2*t**2) for t,v in zip(taus,Ps)]
c3 = neville([(t, v/t**3) for t,v in rem])
check("B.squircle c2 = sqrt2 zeta(2)", c2, mp.sqrt(2)*z2, '1e-10', "NUMERICAL(10d)")
check("B.squircle c3 = sqrt2 zeta(3)", c3, mp.sqrt(2)*z3, '1e-8', "NUMERICAL(8d)")
def P_paper(p):
    eta = 1/p - 1
    f = lambda t: (t*(1-t))**(eta/2)*mp.sqrt(mp.cosh(eta*mp.log(t/(1-t))))
    return 4*mp.sqrt(2)*(1+eta)*mp.quad(f, [0, mp.mpf('0.5'), 1])
h = mp.mpf('0.004')
P3sq = (P_paper(1+2*h) - 2*P_paper(1+h) + 2*P_paper(1-h) - P_paper(1-2*h))/(2*h**3)
check("B.square P3(1) = 6 sqrt2 z3 - 2 sqrt2 pi^2", P3sq, 6*mp.sqrt(2)*z3 - 2*mp.sqrt(2)*pi2, '2e-6', "NUMERICAL(6d)")

# ---------- Task B: C_m law values ----------
def C_m(m):
    s = mp.sin(mp.pi/(2*m)); c = mp.cos(mp.pi/(2*m))
    return 48*m*s**4*c*(3*s**2-1)
check("B.C_2 = 6 sqrt2", C_m(2), 6*mp.sqrt(2), '1e-45', "EXACT")
check("B.C_3 = -9 sqrt3/8", C_m(3), -9*mp.sqrt(3)/8, '1e-45', "EXACT")
check("B.C_6 formula vs paper", C_m(6), mp.mpf('-0.9974443605274988364594'), '1e-12', "paper")

# ---------- Task B: m=3 hexagon via exact-jet engine ----------
exec(open("step_18.py").read().split("print(")[0].replace("mp.mp.dps = 70","mp.mp.dps = 45"))
ln2 = mp.log(2)
P1m3 = P_deriv(3, 1)
check("B.hexagon P1(1) = 2 sqrt3 ln2", P1m3, 2*mp.sqrt(3)*ln2, '1e-35', "NUMERICAL(35d)")
P3m3 = P_deriv(3, 3)
closed = -9*mp.sqrt(3)/8*z3 + mp.sqrt(3)*(3*pi2*ln2/8 + 2*ln2**3 - 3*pi2/4 - 12*ln2**2 + 12*ln2)
check("B.hexagon P3(1) with zeta(3)-coefficient C_3", P3m3, closed, '1e-30', "NUMERICAL(30d)")

# ---------- Task C: the two new families ----------
VS = [(1,1),(1,-1),(-1,1),(-1,-1)]
def hF1(th, tau):
    c, s = mp.cos(th), mp.sin(th); m_ = abs(c)+abs(s)
    L = mp.fsum(mp.exp((v[0]*c+v[1]*s-m_)/tau) for v in VS)
    return m_*mp.exp(tau*mp.log(L))
def hF2(th, tau):
    c, s = mp.cos(th), mp.sin(th)
    return tau*mp.log(mp.fsum(mp.exp((v[0]*c+v[1]*s)/tau) for v in VS))
def P_of(hfun, tau):
    pts = {mp.mpf(0), 2*mp.pi}
    for cth in [mp.mpf(0), mp.pi/2, mp.pi, 3*mp.pi/2]:
        pts.add((cth-20*tau) % (2*mp.pi)); pts.add((cth+20*tau) % (2*mp.pi))
    pts = sorted(pts)
    return mp.fsum(mp.quad(lambda th: hfun(th, tau), [a,b]) for a,b in zip(pts[:-1], pts[1:]) if b-a > 1e-30)
taus = [mp.mpf('0.02')/2**k for k in range(5)]
P1s = [P_of(hF1, t)-8 for t in taus]; P2s = [P_of(hF2, t)-8 for t in taus]
c2a = neville([(t,v/t**2) for t,v in zip(taus,P1s)])
r1 = [(t, v - c2a*t**2) for t,v in zip(taus,P1s)]
c3a = neville([(t,v/t**3) for t,v in r1])
check("C.F1 c2 = pi^2/3", c2a, pi2/3, '1e-10', "NUMERICAL(10d)")
check("C.F1 c3 = 2 zeta(3) (NEW family with zeta(3))", c3a, 2*z3, '1e-8', "NUMERICAL(8d)")
c2b = neville([(t,v/t**2) for t,v in zip(taus,P2s)])
r2 = [(t, v - c2b*t**2) for t,v in zip(taus,P2s)]
c3b = neville([(t,v/t**3) for t,v in r2])
check("C.F2 c3 = 0 (no zeta(3))", c3b, mp.mpf(0), '1e-8', "NUMERICAL")

npass = sum(1 for r in results if r[1])
print()
print("verify_av5b: %d/%d checks PASS in %.1fs" % (npass, len(results), time.time()-t_start))
print("OVERALL:", "PASS" if npass == len(results) else "FAIL")
'''
with open("verify_av5b.py","w") as f: f.write(ver)
print("written OK,", len(ver), "chars")
import py_compile
py_compile.compile("verify_av5b.py", doraise=True)
print("compiles OK")
