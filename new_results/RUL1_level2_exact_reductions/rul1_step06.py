
import sys, time, pickle
sys.dont_write_bytecode = True
import sympy as sp

def log(s):
    s = time.strftime('%H:%M:%S ') + s
    print(s, flush=True)
    open('progress.log','a').write(s+'\n')

t0 = time.time()
log("rul1_step06: TARGET 3 -- regular hexagon (r=1) from thm:hexagon5")

K1G = '<project>/research/08_theory/K1g_hexagon_order5_assembly'
K1F = '<project>/research/08_theory/K1f_hexagon_order5'
bb = pickle.load(open(f'{K1G}/build_blocks.pkl','rb'))
JLOC_n = bb['JLOC_n']
E5 = pickle.load(open(f'{K1F}/E5_polys.pkl','rb'))[5]

A0,A1,B,kap,L2,q0,V02,lA0,lA1,lB = sp.symbols('A0 A1 B kap L2 q0 V02 lA0 lA1 lB')
l2, l3, pi2, z3, z5, sq3 = sp.symbols('l2 l3 pi2 z3 z5 sq3')
z2s, z3s, z5s = sp.symbols('z2 z3 z5')

# paper's regular hexagon: K_p = {sum_j |<n_j,x>|^p <= 1} (w_j = 1), edge length 1/sqrt(3).
# edge V0 = (1/2, 1/(2 sqrt3)) -> V1 = (0, 1/sqrt3): A0=A1=B=1/2, kap=0, L2=1/3, q0=-1/6, V02=1/3
# logs: log(A_j/w_j) = log(1/2) = -log2 (no log3 anywhere)
rd = {A0: sp.Rational(1,2), A1: sp.Rational(1,2), B: sp.Rational(1,2), kap: sp.Integer(0),
      L2: sp.Rational(1,3), q0: sp.Rational(-1,6), V02: sp.Rational(1,3),
      lA0: -l2, lA1: -l2, lB: -l2}

Earc = sp.Integer(0)
for (n,a,b,c), co in E5.items():
    if c != 0: continue
    cov = co.subs(rd)
    if cov == 0: continue
    jv = JLOC_n[(n,a,b)].subs({z2s: pi2/6, z3s: z3, z5s: z5})
    Earc += cov*jv
Earc = sp.expand(120*Earc)
P5 = sp.expand(6*Earc)   # P^(5)(1) = 2 * sum over 3 edges (equal)
log("E_arc and P5 assembled (%.1fs)" % (time.time()-t0))

# check log3-free (logs contain no log3 by construction here)
import collections
cnt = collections.Counter()
for t in sp.Add.make_args(P5):
    p = sp.degree(t, l3) if t.has(l3) else 0
    cnt[p] += 1
print("P5 terms by log3-power:", dict(cnt))
l3free = all(k == 0 for k in cnt)
print("log3-free:", l3free)
assert l3free, "log3 does not cancel in P5"
log("P5 is log3-free (as the paper's closed form)")

# paper's 15-term closed form: (128/sqrt3) P5 = RHS
RHS = (30720*l2 - 4920*pi2 - 61440*l2**2 + 6240*pi2*l2 + 30720*l2**3 + 448*pi2**2
       + 5760*z3*l2 - 1920*pi2*l2**2 - 5120*l2**4 + 11340*z5 + 5580*pi2*z3/6
       - 123*pi2**2*l2 - 1440*z3*l2**2 + 160*pi2*l2**3 + 256*l2**5)
LHS = sp.expand(128/sp.sqrt(3)*P5)
diff = sp.expand(LHS - RHS)
print("(128/sqrt3)P5 - paper's 15-term form =", sp.expand(diff))
assert diff == 0
log("TARGET 3a PASS: (128/sqrt3) P^(5)(1) == paper's 15-term closed form EXACTLY")

# zeta(5), zeta(2)zeta(3) coefficients
P5e = sp.expand(P5)
cz5 = sp.Integer(0)
for t in sp.Add.make_args(P5e):
    if t.has(z5) and not t.has(z3) and not t.has(pi2) and not t.has(l2) and not t.has(l3):
        cz5 += t
cz5 = sp.simplify(cz5/z5)
print("z5 coefficient of P5 =", cz5)
assert sp.simplify(cz5 - 2835*sp.sqrt(3)/32) == 0
cpi2z3 = sp.Integer(0)
for t in sp.Add.make_args(P5e):
    if t.has(z3) and t.has(pi2) and not t.has(z5) and not t.has(l2) and not t.has(l3):
        cpi2z3 += t
cpi2z3 = sp.simplify(cpi2z3/(pi2*z3))
print("pi2*z3 coefficient of P5 =", cpi2z3, "; as z2z3 = pi2 z3/6:", sp.simplify(6*cpi2z3))
assert sp.simplify(6*cpi2z3 - 1395*sp.sqrt(3)/32) == 0
log("TARGET 3b PASS: zeta5 = 2835 sqrt3/32, zeta2zeta3 = 1395 sqrt3/32 exactly")

# numeric check vs the paper's direct value
import mpmath as mp
mp.mp.dps = 60
val = P5.subs({l2: mp.log(2), pi2: mp.pi**2, z3: mp.zeta(3), z5: mp.zeta(5), sq3: mp.sqrt(3)})
val = mp.mpf(sp.N(val, 60))
print("P5 =", mp.nstr(val, 40))
print("paper: 656.4153619676247544766476...")
print("diff =", mp.nstr(abs(val - mp.mpf("656.4153619676247544766476")), 3))
log("TARGET 3 numeric agreement printed; step06 done %.1fs" % (time.time()-t0))
