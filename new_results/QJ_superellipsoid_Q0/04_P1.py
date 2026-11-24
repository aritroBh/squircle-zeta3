# 04_P1.py -- P1 = int_0^1 [2 log(1+q) - (1+q) log(2q)]/(1-q^2) ds,  q = sqrt(s^2+(1-s)^2).
# Euler substitution q = (v^2+1)/(2 sqrt2 v), v in [1, a]; partial fractions (exact residues);
# six primitive types, each verified by exact differentiation under its stated domain hypothesis;
# the hypotheses are checked exactly for every term; complex letter log(v^2+1) handled by Re Li2 along
# segments that avoid the cut [1,oo); its three boundary values by reflection/inversion/Landen + Lemma A.
import sympy as sp, mpmath as mp, json
from common import lam, l3, L, al, PI, Li2, li2, logS, num, show, Q2, A, B, SQ2, q2
mp.mp.dps = 90
ref = json.load(open("reference_100.json"))
def ok(e):
    if sp.simplify(e) == 0: return True
    return sp.radsimp(sp.simplify(sp.expand(sp.sympify(e)))) == 0
v, s, y, d, D, z, X = sp.symbols('v s y d D z X', positive=True)
r2 = sp.sqrt(2); a_ = 1 + r2; b_ = r2 - 1

print("== Step 1: Euler substitution")
qv = (v**2 + 1)/(2*r2*v)
assert ok(2*qv**2 - 1 - ((v**2 - 1)/(2*v))**2)
assert ok(qv*sp.diff(qv, v)/((v**2 - 1)/(2*v)) - (v**2 + 1)/(4*v**2))
assert ok(qv.subs(v, 1) - 1/r2) and ok(qv.subs(v, a_) - 1)
assert ok(sp.diff(qv, v) - (v**2 - 1)/(2*r2*v**2))          # q increasing on [1,a]
print("  q=(v^2+1)/(2sqrt2 v): sqrt(2q^2-1)=(v^2-1)/(2v) (v>=1), q dq/sqrt(2q^2-1)=(v^2+1)/(4v^2) dv, q(1)=1/sqrt2, q(a)=1  (exact)")
print("  => P1 = 2 int_{1/sqrt2}^1 f(q) q dq/sqrt(2q^2-1) = (1/2) int_1^a f(q(v)) (v^2+1)/v^2 dv")
assert ok(1 + qv - (v + a_)*(v + b_)/(2*r2*v)) and ok(1 - qv + (v - a_)*(v - b_)/(2*r2*v)) and ok(2*qv - (v**2 + 1)/(r2*v))
print("  1+q=(v+a)(v+b)/(2sqrt2 v), 1-q=-(v-a)(v-b)/(2sqrt2 v), 2q=(v^2+1)/(sqrt2 v)  (exact)")
fq = lambda qq: (2*mp.log(1 + qq) - (1 + qq)*mp.log(2*qq))/(1 - qq*qq)
numv = mp.quad(lambda t: fq((t*t + 1)/(2*mp.sqrt(2)*t))*(t*t + 1)/(2*t*t), [1, 2, 1 + mp.sqrt(2)])
show("  P1 (v-form, quadrature) vs reference", numv, mp.mpf(ref["P1"]), 80)

# integrand  (v^2+1)/(2v^2) f(q)  =  sum_letters C_l(v) log l(v) + C_lam(v) * lambda
# 2 log(1+q) = 2log(v+a) + 2log(v+b) - 2 log v - 3 lambda ;  log(2q) = log(v^2+1) - log v - lambda/2
pref = (v**2 + 1)/(2*v**2*(1 - qv**2))
C = {'v+a': 2*pref, 'v+b': 2*pref, 'v': -2*pref + pref*(1 + qv), 'v2+1': -pref*(1 + qv),
     'lam': -3*pref + pref*(1 + qv)/2}
lv = sp.Function('lg')
recon = sum(C[k]*lv(k) for k in C if k != 'lam') + C['lam']*lam
lam_expr = {'v+a': sp.log(v + a_), 'v+b': sp.log(v + b_), 'v': sp.log(v), 'v2+1': sp.log(v**2 + 1)}
orig = (v**2 + 1)/(2*v**2)*(2*(sp.log(v + a_) + sp.log(v + b_) - sp.log(v) - sp.Rational(3, 2)*sp.log(2))
                             - (1 + qv)*(sp.log(v**2 + 1) - sp.log(v) - sp.log(2)/2))/(1 - qv**2)
chk = recon.subs({lv(k): lam_expr[k] for k in lam_expr}).subs(lam, sp.log(2)) - orig
assert ok(chk)
# and orig equals (v^2+1)/(2v^2) f(q(v)) after expanding logs of positive quantities:
fv = (v**2 + 1)/(2*v**2)*(2*sp.log(1 + qv) - (1 + qv)*sp.log(2*qv))/(1 - qv**2)
for vv in [sp.Rational(11, 10), sp.Rational(3, 2), sp.Rational(23, 10)]:
    assert abs(sp.N((orig - fv).subs(v, vv), 60)) < 1e-55
print("  letters: log(v+a), log(v+b), log v, log(v^2+1), lambda; coefficient functions C_l (exact)")

print("== Step 2: partial fractions (exact residues at the simple poles 0, +-a, +-b)")
poles = {'0': Q2(0), 'a': A, '-a': -A, 'b': B, '-b': -B}
PF = {}
for k, Ck in C.items():
    Ck = sp.cancel(sp.together(Ck))
    res = {}
    for pn, pv in poles.items():
        rr = sp.radsimp(sp.simplify(sp.limit((v - pv.sym())*Ck, v, pv.sym())))
        if rr != 0: res[pn] = q2(rr)
    rec = sum(res[pn].sym()/(v - poles[pn].sym()) for pn in res)
    assert ok(Ck - rec), k
    PF[k] = res
    print("   %-5s" % k, {pn: res[pn] for pn in res})

# residue cancellation at v = a:  sum_l r_{l,a} log l(a) + r_{lam,a} lambda = 0
la = {'v+a': logS(2*A), 'v+b': logS(A + B), 'v': 0, 'v2+1': logS(A*A + 1)}
tot = sum(PF[k].get('a', Q2(0)).sym()*la[k] for k in la) + PF['lam'].get('a', Q2(0)).sym()*lam
assert sp.expand(tot) == 0
print("  residues at v=a cancel exactly: sum_l r_(l,a) log l(a) + r_(lam,a) lambda = 0")

print("== Step 3: primitive types (each verified by exact differentiation under its hypothesis)")
# T1: int [log(v+p)-log(a+p)]/(v-a) dv = -Li2((a-v)/(a+p))            (v+p>0, a+p>0, v<=a)
#     with X = (a-v)/(a+p) in [0,1):  d/dv = -(-log(1-X)/X)(-1/(a+p)) ... check in X:
assert ok(sp.expand_func(sp.diff(-sp.polylog(2, X), X)) - sp.log(1 - X)/X)
print("  T1: -Li2(X), X=(a-v)/(a+p): d/dX = log(1-X)/X, dX/dv=-1/(a+p), v-a=-(a+p)X, log(1-X)=log(v+p)-log(a+p)  (exact)")
# T2: d=c+p>0, y=v-c>0:  int log(v+p)/(v-c) dv = log d log y - Li2(-y/d)
F2 = sp.log(d)*sp.log(y) - sp.polylog(2, -y/d)
assert ok(sp.expand_log(sp.expand_func(sp.diff(F2, y)) - sp.log(y + d)/y, force=True))
print("  T2: d>0,y>0: d/dy[log d log y - Li2(-y/d)] = log(y+d)/y  (exact; logs of positive numbers)")
# T3: d=c+p<0, D=-d, y=v-c>D>0 (y = D + z):  int = (1/2)log^2 y + Li2(D/y)
F3 = sp.log(y)**2/2 + sp.polylog(2, D/y)
e3 = (sp.expand_func(sp.diff(F3, y)) - sp.log(y - D)/y).subs(y, D + z)
assert ok(sp.expand_log(sp.together(e3), force=True)) or ok(sp.expand_log(sp.simplify(e3), force=True))
print("  T3: d<0,y>|d|: d/dy[(1/2)log^2 y + Li2(|d|/y)] = log(y-|d|)/y  (exact)")
# T4: d=0: int log(y)/y = (1/2) log^2 y   (trivial)
# T5 (complex letter, pole c=0): d/dv[-(1/2)Li2(-v^2)] = log(1+v^2)/v
assert ok(sp.expand_func(sp.diff(-sp.polylog(2, -v**2)/2, v)) - sp.log(1 + v**2)/v)
print("  T5: d/dv[-(1/2)Li2(-v^2)] = log(1+v^2)/v  (exact)")
# T6 (complex letter, pole c real): Z(v) = (c-v)/(c+i):  d/dv Li2(Z) = -log(1-Z) Z'/Z, Z'/Z = 1/(v-c) (real),
#     |1-Z|^2 = (v^2+1)/(c^2+1).  Z(v) is real only at v=c (Z=0), so Z(v) never meets [1,oo): Li2(Z(v)) analytic.
#     => d/dv[-2 Re Li2(Z(v))] = [log(v^2+1) - log(c^2+1)]/(v-c)
cc = sp.symbols('c', real=True); vr = sp.symbols('v', real=True)
Zc = (cc - vr)/(cc + sp.I)
assert ok(sp.diff(Zc, vr)/Zc - 1/(vr - cc))
assert ok(sp.expand((1 - Zc)*(1 - sp.conjugate(Zc))) - (vr**2 + 1)/(cc**2 + 1))
assert ok(sp.expand_func(sp.diff(sp.polylog(2, X), X)) + sp.log(1 - X)/X)
print("  T6: Z=(c-v)/(c+i): Z'/Z = 1/(v-c), |1-Z|^2=(v^2+1)/(c^2+1)  (exact) => d/dv[-2 Re Li2(Z)] = [log(v^2+1)-log(c^2+1)]/(v-c)")

print("== Step 4: assemble P1 = sum over (letter, pole) of primitive differences on [1, a]")
ReZ = {}                       # symbolic placeholders for Re Li2 at the three complex points
Zs = sp.symbols('ReLi2_Z1 ReLi2_Z2 ReLi2_Z3')
v0, v1 = Q2(1), A
P1 = 0
contrib = []
def add(lbl, e):
    global P1
    P1 += e; contrib.append((lbl, sp.expand(e)))
pofl = {'v+a': A, 'v+b': B, 'v': Q2(0)}
for k, res in PF.items():
    for pn, r in res.items():
        c = poles[pn]; R = r.sym()
        if k == 'lam':
            if pn == 'a': continue                       # cancelled (Step 2)
            assert c < v0
            add("lambda*%s/(v-%r) [log]" % (r, c), R*lam*(logS(v1 - c) - logS(v0 - c)))
            continue
        if k == 'v2+1':
            if pn == '0':
                add("log(v^2+1)*%s/v [T5]" % r, R*(-li2(-(v1*v1))/2 + li2(-(v0*v0))/2))
            elif pn == 'a':
                # [log(v^2+1)-log(a^2+1)]/(v-a): -2Re Li2((a-v)/(a+i)); at v=a: 0; at v=1: Z1=(a-1)/(a+i)
                add("log(v^2+1)*%s/(v-a) [T6]" % r, R*(0 + 2*Zs[0]))
            elif pn == 'b':
                # log(v^2+1)/(v-b) = log(b^2+1)/(v-b) + [..]: log(b^2+1) log(v-b) - 2 Re Li2((b-v)/(b+i))
                assert c < v0
                add("log(v^2+1)*%s/(v-b) [T6]" % r, R*(logS(B*B + 1)*(logS(v1 - c) - logS(v0 - c)) - 2*Zs[1] + 2*Zs[2]))
            else: raise ValueError
            continue
        p = pofl[k]
        assert (v0 + p) > 0 and (v1 + p) > 0
        if pn == 'a':
            assert A + p > 0
            Xv0 = (A - v0)/(A + p); assert Xv0 >= 0 and Xv0 < 1
            add("log(%s)*%s/(v-a) [T1]" % (k, r), R*(0 - (-li2(Xv0))))   # T1 (constant part cancelled in Step 2)
            continue
        assert c < v0                                       # y = v - c > 0 on [1,a]
        dd = c + p
        typ = "T2" if dd > 0 else ("T3" if dd < 0 else "T4")
        if dd > 0:   # T2
            F = lambda vv: logS(dd)*logS(vv - c) - li2(-((vv - c)/dd))
        elif dd < 0: # T3
            assert (v0 - c) > -dd
            F = lambda vv: logS(vv - c)**2/2 + li2((-dd)/(vv - c))
        else:        # T4
            F = lambda vv: logS(vv - c)**2/2
        add("log(%s)*%s/(v-%r) [%s]" % (k, r, c, typ), R*(F(v1) - F(v0)))
P1 = sp.expand(P1)
print("  term-by-term contributions (letter * residue/(v-c) [primitive type]):")
for lbl, e in contrib: print("    %-34s -> %s" % (lbl, e))

print("== Step 5: the three complex boundary values (conjugates of z1, z2, z3 of the notes; Re is conjugation-invariant)")
I = sp.I
Z1 = (a_ - 1)/(a_ + I); Z2 = (b_ - a_)/(b_ + I); Z3 = (b_ - 1)/(b_ + I)
assert ok(sp.radsimp(Z1) - (1 - I*b_)/2) and ok(sp.radsimp(Z2) - (-1 + I*a_)/r2) and ok(sp.radsimp(Z3) - (-b_ + I)/2)
print("  Z1=(a-1)/(a+i)=(1-ib)/2, Z2=(b-a)/(b+i)=(-1+ia)/sqrt2, Z3=(b-1)/(b+i)=(-b+i)/2  (exact)")
# exact angles: tan(pi/8) = b, tan(3pi/8) = a
assert ok(sp.tan(PI/8) - b_) and ok(sp.tan(3*PI/8) - a_)
LiN = lambda zz: mp.polylog(2, zz)
bN, aN = mp.sqrt(2) - 1, mp.sqrt(2) + 1
LemA = li2(B)/2 + PI**2/192                                  # Re Li2(b e^{+-i pi/4}), Lemma A
# Z1: Re Z1 = 1/2, 1-Z1 = conj(Z1); reflection Li2(z)+Li2(1-z) = pi^2/6 - log z log(1-z) (z notin (-oo,0]U[1,oo))
#     => 2 Re Li2(Z1) = pi^2/6 - |log Z1|^2,  |Z1|^2 = b/sqrt2, arg Z1 = -pi/8
assert ok(sp.expand(Z1*sp.conjugate(Z1)) - b_/r2) and ok(1 - Z1 - sp.conjugate(Z1))
lZ1 = logS(B/SQ2)/2
vZ1 = sp.expand(PI**2/12 - (lZ1**2 + (PI/8)**2)/2)
z1N = (1 - 1j*bN)/2
show("  Re Li2(Z1) [reflection]", mp.re(LiN(z1N)), num(vZ1), 85)
# Z2: 1-Z2 = a e^{-i pi/4}; reflection then inversion Li2(w)+Li2(1/w) = -pi^2/6 - log^2(-w)/2 (w notin [0,oo)),
#     w = a e^{-i pi/4}, 1/w = b e^{i pi/4}, -w = a e^{3 i pi/4};  |Z2|^2 = sqrt2 a, arg Z2 = pi - 3pi/8 = 5pi/8
assert ok(sp.radsimp(1 - Z2) - a_*(1 - I)/r2) and ok(sp.expand(Z2*sp.conjugate(Z2)) - r2*a_)
ReLi_w = -PI**2/6 - (L**2 - (3*PI/4)**2)/2 - LemA
vZ2 = sp.expand(PI**2/6 - (logS(SQ2*A)/2*L - (5*PI/8)*(-PI/4)) - ReLi_w)
z2N = (-1 + 1j*aN)/mp.sqrt(2)
show("  Re Li2(a e^{-i pi/4}) [inversion + Lemma A]", mp.re(LiN(aN*mp.expjpi(-mp.mpf(1)/4))), num(ReLi_w), 85)
show("  Re Li2(Z2) [reflection]", mp.re(LiN(z2N)), num(vZ2), 85)
# Z3: Landen Li2(z) + Li2(z/(z-1)) = -log^2(1-z)/2 (z notin [1,oo)); Z3/(Z3-1) = b e^{-i pi/4};
#     |1-Z3|^2 = a/sqrt2, arg(1-Z3) = -pi/8
assert ok(sp.radsimp(Z3/(Z3 - 1)) - b_*(1 - I)/r2)
assert ok(sp.expand((1 - Z3)*sp.conjugate(1 - Z3)) - a_/r2)
vZ3 = sp.expand(-((logS(A/SQ2)/2)**2 - (PI/8)**2)/2 - LemA)
z3N = (-bN + 1j)/2
show("  Re Li2(Z3) [Landen + Lemma A]", mp.re(LiN(z3N)), num(vZ3), 85)
for nm, zz in (("Z1", z1N), ("Z2", z2N), ("Z3", z3N)):
    print("     %s = %s, arg = %s*pi" % (nm, mp.nstr(zz, 8), mp.nstr(mp.arg(zz)/mp.pi, 8)))
print("     arg(1-Z3)/pi =", mp.nstr(mp.arg(1 - z3N)/mp.pi, 10), " arg Z2/pi =", mp.nstr(mp.arg(z2N)/mp.pi, 10))
P1 = sp.expand(P1.subs({Zs[0]: vZ1, Zs[1]: vZ2, Zs[2]: vZ3}))
print("== Result (unreduced):")
print("  Li2 arguments:", sorted(set(str(f.args[0]) for f in P1.atoms(sp.Function))))
print("  P1 =", P1)
show("P1 exact (unreduced) vs 100-digit reference", num(P1), mp.mpf(ref["P1"]), 90)
json.dump({"P1": sp.srepr(P1)}, open("P1_unreduced.json", "w"))
