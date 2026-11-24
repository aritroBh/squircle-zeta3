# Independent referee checks for thm:M1 / lem:dilog / app:M1 (own code).
import mpmath as mp, sympy as sp
mp.mp.dps = 110
s2 = mp.sqrt(2); a = 1+s2; b = s2-1; L = mp.log(a); lam = mp.log(2)
Li2 = lambda z: mp.re(mp.polylog(2, z))   # real part convention for z>1
pi2 = mp.pi**2

def closed():
    return (19*s2/48-mp.mpf(1)/24)*pi2 + lam**2/2 + 2*s2*lam*L - (1+5*s2/2)*L**2 - 5*s2*Li2(b)

C = closed()
# Route A: eq (M1t), t-form
q = lambda t: mp.sqrt(t**2+(1-t)**2)
fA = lambda t: mp.log(t)/t*((1-t)/q(t)-1)
A = mp.quad(fA, [0, mp.mpf(1)/4, mp.mpf(1)/2, mp.mpf(3)/4, 1])
# Route A': unsymmetrized derivative of Phi at eps=0, from scratch:
# dPhi/deps = log t * a + log(1-t)*b - (log t a^2 + log(1-t) b^2)/sqrt(a^2+b^2), a=1/t,b=1/(1-t)
def dPhi(t):
    x, y = mp.log(t), mp.log(1-t); A_, B_ = 1/t, 1/(1-t)
    return x*A_ + y*B_ - (x*A_**2 + y*B_**2)/mp.sqrt(A_**2+B_**2)
Fp = mp.quad(dPhi, [0, mp.mpf(1)/4, mp.mpf(1)/2, mp.mpf(3)/4, 1])
Ap = -Fp/2
# Route C: v-form from sec:method
cst = 1-1/s2
fC = lambda v: (mp.log(v-b)+mp.log(v+a)-mp.log(v)-2*lam)*(cst/v-2/(v+a))
Cv = mp.quad(fC, [b, 1, a])
print("closed form      M1 =", mp.nstr(C, 100))
print("route A (M1t)    M1 =", mp.nstr(A, 100))
print("route A' (-F'/2) M1 =", mp.nstr(Ap, 100))
print("route C (v-form) M1 =", mp.nstr(Cv, 100))
for nm, val in [("A", A), ("A'", Ap), ("C", Cv)]:
    print(f"|closed - {nm}| =", mp.nstr(abs(C-val), 5))
print("kappa = 8 M1 =", mp.nstr(8*C, 40), " printed: 2.02095276965635835695...")
print("Ahmed 10 digits 2.0209527696 ; 8*M0 =", mp.nstr(8*(lam-s2*L),20), " -4K =", mp.nstr(-4*(2*s2*L-2*lam),20))

# sympy: v-substitution and rational integrand
v = sp.symbols('v', positive=True); S2 = sp.sqrt(2); A_s = 1+S2; B_s = S2-1
t = (v**2+2*v-1)/(4*v)
print("t factor:", sp.simplify(t-(v-B_s)*(v+A_s)/(4*v)), " 1-t:", sp.simplify(1-t-(A_s-v)*(v+B_s)/(4*v)))
qv = (v**2+1)/(2*S2*v)
print("q^2 check:", sp.simplify(sp.expand(t**2+(1-t)**2-qv**2)), " dt:", sp.simplify(sp.diff(t,v)-(v**2+1)/(4*v**2)))
print("t(b)=", sp.nsimplify(sp.simplify(t.subs(v,B_s))), " t(a)=", sp.simplify(t.subs(v,A_s)), " q>0 on [b,a]: yes (v>0)")
rat = ((1-t)/qv - 1)/t*sp.diff(t,v)
print("rational kernel - (c*/v - 2/(v+a)) =", sp.simplify(sp.radsimp(rat-((1-1/S2)/v-2/(v+A_s)))))

# I(x,y) table
def I(x, y): return mp.quad(lambda v: mp.log(abs(v-x))/(v-y), [b, 1, a])
tab = {
 "I(b,0)":  (I(b,0),  Li2(b**2)-pi2/6),
 "I(-a,0)": (I(-a,0), Li2(-b**2)+pi2/12+2*L**2),
 "I(0,0)":  (I(0,0),  mp.mpf(0)),
 "I(-a,-a)":(I(-a,-a),-mp.mpf(5)/8*lam**2+lam*L+L**2/2),
 "I(b,-a)": (I(b,-a), Li2(2-s2)-pi2/6-mp.mpf(5)/8*lam**2+lam*L+L**2/2),
 "I(0,-a)": (I(0,-a), -Li2((2+s2)/4)+pi2/12-mp.mpf(9)/8*lam**2+lam*L+L**2/2)}
for k,(n,f) in tab.items(): print(k, "diff =", mp.nstr(abs(n-f),5))
print("int dv/v - 2L:", mp.nstr(mp.quad(lambda v:1/v,[b,a])-2*L,5), " int dv/(v+a) - (L-lam/2):", mp.nstr(mp.quad(lambda v:1/(v+a),[b,a])-(L-lam/2),5))
# w-endpoints: show which I's need dilog at w>1 (branch cut / inversion for z>1)
for (x,y,nm) in [(b,0,"I(b,0)"),(-a,0,"I(-a,0)"),(b,-a,"I(b,-a)"),(0,-a,"I(0,-a)")]:
    print(nm, "w0 =", mp.nstr((b-y)/(x-y),8), " w1 =", mp.nstr((a-y)/(x-y),8))
# assembly
Ivals = {k: v_[0] for k,v_ in tab.items()}
asm = cst*(Ivals["I(b,0)"]+Ivals["I(-a,0)"]-Ivals["I(0,0)"]) - 2*(Ivals["I(b,-a)"]+Ivals["I(-a,-a)"]-Ivals["I(0,-a)"]) - 2*lam**2 + 2*s2*lam*L
print("assembly - closed:", mp.nstr(asm-C,5))
br1 = Ivals["I(b,0)"]+Ivals["I(-a,0)"]-Ivals["I(0,0)"]
print("c* bracket - (10Li2(b)-19pi^2/24+5L^2):", mp.nstr(br1-(10*Li2(b)-19*pi2/24+5*L**2),5))
rest = -2*(Ivals["I(b,-a)"]+Ivals["I(-a,-a)"]-Ivals["I(0,-a)"]) - 2*lam**2 + 2*s2*lam*L
print("rest - (3pi^2/4+lam^2/2+2sqrt2 lam L-6L^2-10Li2(b)):", mp.nstr(rest-(3*pi2/4+lam**2/2+2*s2*lam*L-6*L**2-10*Li2(b)),5))

# Lemma dilog lines
lines = [
 ("Li2(b)-Li2(-b)", Li2(b)-Li2(-b), pi2/8-L**2/2),
 ("Li2(b^2)", Li2(b**2), 4*Li2(b)-pi2/4+L**2),
 ("Li2(2-sqrt2)", Li2(2-s2), -Li2(b)+pi2/6+lam*L/2-L**2),
 ("Li2(-sqrt2)", Li2(-s2), Li2(b)-pi2/6-lam*L/2+L**2/2),
 ("Li2(-b^2)", Li2(-b**2), 6*Li2(b)-11*pi2/24+2*L**2),
 ("Li2((2+sqrt2)/4)", Li2((2+s2)/4), 6*Li2(b)-7*pi2/24-mp.mpf(9)/8*lam**2-mp.mpf(3)/2*lam*L+mp.mpf(7)/2*L**2)]
for n,l,r in lines: print("lemma", n, "diff =", mp.nstr(abs(l-r),5))
# identities used, at the stated points, with domain checks
x=-s2; y=-b
print("five-term domain x<1,y<1,x+y<1:", x<1, y<1, x+y<1, " args:", [mp.nstr(z,12) for z in (x/(1-y), y/(1-x), x*y/((1-x)*(1-y)))])
ft = Li2(x/(1-y))+Li2(y/(1-x))-Li2(x*y/((1-x)*(1-y))) - (Li2(x)+Li2(y)+mp.log(1-x)*mp.log(1-y))
print("five-term residual:", mp.nstr(ft,5))
chi = lambda u: (Li2(u)-Li2(-u))/2
print("chi2 Landen at b residual:", mp.nstr(chi((1-b)/(1+b))+chi(b)-pi2/8-mp.log(b)*mp.log((1+b)/(1-b))/2,5), " fixed point:", mp.nstr((1-b)/(1+b)-b,5))
z=(2+s2)/4; print("Landen z=(2+s2)/4: z/(z-1)+a^2 =", mp.nstr(z/(z-1)+a**2,5), " residual:", mp.nstr(Li2(z)+Li2(z/(z-1))+mp.log(1-z)**2/2,5))
z=-s2; print("Landen z=-sqrt2: z/(z-1)=", mp.nstr(z/(z-1),12), " residual:", mp.nstr(Li2(z)+Li2(z/(z-1))+mp.log(1-z)**2/2,5))
print("inversion z=-a^2 residual:", mp.nstr(Li2(-a**2)+Li2(-b**2)+pi2/6+mp.log(a**2)**2/2,5))
# inversion for z>1 (needed for I(b,0), I(b,-a), I(0,-a)):
for zz in (a**2, a/s2, 2*s2*b):
    print("inversion z>1 at", mp.nstr(zz,8), ": ReLi2(z)+Li2(1/z)-pi^2/3+log^2 z/2 =", mp.nstr(Li2(zz)+Li2(1/zz)-pi2/3+mp.log(zz)**2/2,5))
