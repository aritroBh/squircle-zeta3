# step_10.py - D13h FINAL ANALYTIC CELL: one C for ALL N >= 1e10 (y in [0, y_J])
# Every piece = y^3 * (regular factor); regular factors evaluated on the box [0, y_J]
# with NO division by y. All arb enclosures. Derivation in comments.
from flint import arb, acb, ctx
import time
ctx.dps = 70
pi = arb.pi(); rt2 = arb(2).sqrt(); sqpi = pi.sqrt()
a = pi/(2*rt2)
def iv(lo, hi):
    lo = arb(lo).lower(); hi = arb(hi).upper()
    return arb((lo+hi)/2, (hi-lo)/2)
def up(x): return arb(x.upper())
def lo(x): return arb(x.lower())

t0 = time.time()
# ---------- box ----------
N_J = arb(10)**10
y_J = (pi/(2*N_J)).sqrt()              # y_J = sqrt(pi/2 * 1e-10)
Y  = iv(0, y_J)                        # y box
SG = iv(0, 1)                          # sigma box (nsub=1; rho-1 = y^2 sig <= 1.6e-10)
R  = iv(1, 1 + y_J**2)                 # rho box
rU = up(R)                             # rho_hat
Tval = up(R*y_J)                       # T = rho_hat*y_J
T  = iv(0, Tval)                       # t = rho*y box
c_w = arb("1.0")
L_J = (-y_J.log())                     # log(1/y_J) = 11.7356...
al_ = (pi/2).log() + c_w               # a_ = log(pi/2)+c_w ; w^2 = al_ + 2L
w_J = (al_ + 2*L_J).sqrt()             # w at y_J (min of w on box)
assert (al_ + 2*L_J) > 6, "monotonicity windows for y*w^q fail"
sqyJ0 = y_J.sqrt()   # sqrt(y_J) point value (box sqrt of [0,y_J] is nan in arb)

# ---------- regular factor box-evals ----------
def numR_f(Tb): return Tb**3 + 2*rt2*Tb**2 + 4*Tb + rt2
def D2_f(Tb):  return Tb**2 + 2*rt2*Tb + 4
def den1_f(Yb, Rb): return 1 + rt2*Rb*Yb + Rb**2*Yb**2
# c2reg = (pi rho^3 /4) numR / (1+rt2 t + t^2)^2   [c2 = y*c2reg, curvature at saddle eta=1]
c2reg = (pi*R**3/4)*numR_f(T)/(1+rt2*T+T**2)**2
c2reg_lo, c2reg_up = lo(c2reg), up(c2reg)
print(f"c2reg  = {c2reg}  (limit pi/(2sqrt2)={float(pi/(2*rt2)):.6f})")

# Omn, omreg: om = y*omreg, omreg = (pi/(2 rt2)) Omn/den1
Omn = rt2*R + R**2*Y - 2*Y*SG - Y**3*SG**2
assert lo(Omn) > 0, "Omn lost sign"
omreg = (pi/(2*rt2))*Omn/den1_f(Y, R)
omreg_up = up(omreg)
print(f"omreg  = {omreg}  (limit pi/2)")

# mSreg: mS = y*mSreg, mSreg = (pi rho^3/4) h3, h3 <= rt2/3
mSreg_up = up(pi*R**3/4)*rt2/3
eRe0_up = (y_J*mSreg_up).exp()         # e^{RePhi0+th} <= this
print(f"mSreg_up = {mSreg_up}, eRe0_up = {eRe0_up}")

# ---- g0 Binet bound: |g0| <= t^{5/2} (sqpi/4) 2^{3/2} e^{ReE}, ReE <= (t/24)(-rt2 + 459 t^2/(1-3T)) <= 0
assert (arb(459)*Tval**2/(1-3*Tval)) < rt2, "ReE<=0 fails"
G0reg = sqpi/4 * 2**arb("1.5")
P0reg = G0reg * eRe0_up                # P0t = t^{5/2} P0reg  (P0t = |g0| e^{RePhi0+th})
print(f"G0reg = {G0reg}, P0reg = {P0reg}")

# ---------- delta / zmax ----------
# d = w/sqrt(2 c2) = w / sqrt(2 y c2reg)  [assert d < x0/3 i.e. y w^2 < 2 c2reg/(9 rho^2)]
# y w^2 = y(al_+2L) is increasing on (0,y_J] since al_+2L_J > 2 -> bound by endpoint:
yw2_up = y_J*(al_ + 2*L_J)
assert yw2_up < 2*c2reg_lo/(9*rU**2), "d < x0/3 fails"
# delta = d*t = w rho sqrt(y)/sqrt(2 c2reg);  delta^2 = w^2 rho^2 y/(2 c2reg) <= y_J(al_+2L_J) rU^2/(2 c2reg_lo)
delta_max = (y_J*(al_+2*L_J)).sqrt()*rU/(2*c2reg_lo).sqrt()
zmax_lo = w_J/(2*c2reg_up).sqrt()       # zmax = w rho/sqrt(2c2reg) >= w_J/sqrt(2 c2reg_up) (rho>=1)
zmax_up_w3 = w_J**3                     # for joint bounds w^3 sqrt(y) <= w_J^3 sqrt(y_J)
print(f"delta_max = {float(delta_max):.6f}, zmax_lo = {float(zmax_lo):.6f}, w_J = {float(w_J):.6f}")

# ---------- curvature curv2reg(eta) = (pi rho^2/4) f/g and its eta-Lipschitz ----------
def cf(Eb, Yb, Rb):  # num(eta)
    return rt2*Eb**3*Rb + 4*Eb**2*Rb**2*Yb + 2*rt2*Eb*Rb**3*Yb**2 + Rb**4*Yb**3
def cg(Eb, Yb, Rb):  # den(eta) = eta^2 (eta^2+rt2 eta rho y + rho^2 y^2)^2
    return Eb**2*(Eb**2 + rt2*Eb*Rb*Yb + Rb**2*Yb**2)**2
def curv2reg(Eb): return (pi*R**2/4)*cf(Eb, Y, R)/cg(Eb, Y, R)
# eta box for the reach: [1-delta_max, 1+delta_max]
ERA = iv(1-delta_max, 1+delta_max)
# Lip_eta: |d/deta curv2reg| = (pi rho^2/4)|f' g - f g'|/g^2
def cfp(Eb, Yb, Rb): return 3*rt2*Eb**2*Rb + 8*Eb*Rb**2*Yb + 2*rt2*Rb**3*Yb**2
def cgp(Eb, Yb, Rb):
    h = Eb**2 + rt2*Eb*Rb*Yb + Rb**2*Yb**2
    return 2*Eb*h*(3*Eb**2 + 2*rt2*Eb*Rb*Yb + Rb**2*Yb**2)
fB = cf(ERA, Y, R); gB = cg(ERA, Y, R); fpB = cfp(ERA, Y, R); gpB = cgp(ERA, Y, R)
diffB = fpB*gB - fB*gpB            # = g^2 d/deta(f/g), NEGATIVE (curv2reg decreasing in eta)
Lip_eta = up(pi*rU**2/4)*max(abs(lo(diffB)), abs(up(diffB)))/lo(gB)**2
print(f"Lip_eta(curv2reg) = {Lip_eta}  (3*limit {(float(3*pi/(2*rt2))):.4f}-ish)")

# ---------- Gh = rho(eta+rt2 t)/(rt2 eta (eta^2 + rt2 eta t + t^2)); -Gh' lower on each side ----------
def Gh(Eb, Tb, Rb): return Rb*(Eb+rt2*Tb)/(rt2*Eb*(Eb**2+rt2*Eb*Tb+Tb**2))
def mGhp_num(Eb, Tb): return 2*Eb**3 + 4*rt2*Eb**2*Tb + 4*Eb*Tb**2 + rt2*Tb**3   # -Gh' = (rho/rt2) num/g^2
def Gh_g(Eb, Tb): return Eb*(Eb**2+rt2*Eb*Tb+Tb**2)
ERa = iv(1, 1+delta_max); ELa = iv(1-delta_max, 1)
g1R = (1/rt2)*lo(mGhp_num(ERa, T))/up(Gh_g(ERa, T))**2 * lo(R)
g1L = (1/rt2)*lo(mGhp_num(ELa, T))/up(Gh_g(ELa, T))**2 * lo(R)
print(f"g1R_lo = {g1R}, g1L_lo = {g1L}  (limit sqrt2={float(rt2):.4f})")

# ---------- lamR, lamL regular: lam >= sqrt(y) w_J lamreg ----------
Q1 = R**2 - SG*(2+Y**2*SG)*(1+rt2*T)
Q1abs = max(abs(lo(Q1)), abs(up(Q1)))
yJ32 = y_J**arb("1.5")
lamRreg = (pi*lo(R)**2/2)*g1R/(2*c2reg_up).sqrt() - a*yJ32*Q1abs/w_J
lamLreg = (pi*lo(R)**2/2)*g1L/(2*c2reg_up).sqrt() - a*yJ32*Q1abs/w_J
assert lamRreg > 0 and lamLreg > 0, "lamreg lost sign"
print(f"lamRreg = {lamRreg}, lamLreg = {lamLreg}")

# tail e^{(c2-c2R)d^2} corrections
E1R = (Lip_eta*sqyJ0*w_J**3*rU/(2*c2reg_lo)**arb("1.5")).exp()
E1L = E1R
print(f"E1R = {E1R}")

# ---------- Ph3reg, Ph4reg, G1reg, G2reg (eta-boxed) ----------
def Ph3reg(Eb): return pi*R**4*(R**2*Y**2 + 3*Eb*R*Y + 3*Eb**2)/(Eb**3*(Eb**2+rt2*Eb*R*Y+R**2*Y**2)**arb("1.5"))
def Ph4reg(Eb): return 3*pi*R**5*(R**3*Y**3 + 4*Eb*R**2*Y**2 + 6*Eb**2*R*Y + 4*Eb**3)/(Eb**4*(Eb**2+rt2*Eb*R*Y+R**2*Y**2)**2)
def G1reg(Eb): return 6*R**2*Y/Eb**2 + (arb("1.5")*rt2+1)*R/Eb
def G2reg(Eb): return 12*rt2*R**3*Y/Eb**3 + 4*R**2/Eb**2

# ---------- W1 ----------
G25 = arb("2.5").gamma(); G15 = arb("1.5").gamma()
m4 = G25/c2reg_lo**arb("2.5"); m2 = G15/c2reg_lo**arb("1.5")
W1reg = rU**arb("2.5")*P0reg*omreg_up*(up(Ph3reg(iv(1,1)))*m4/6 + up(G1reg(iv(1,1)))*m2)
print(f"W1reg = {W1reg}  -> C-part {float((pi/2)**arb('1.5')*arb(W1reg)):.4f}")

# ---------- W2, W3 : z-segmentation on [0, zstar], z = u rho sqrt(y) ----------
J = 192; zstar = arb("3.0")
assert zstar <= zmax_lo, "zstar beyond zmax_lo"
dz = zstar/J
K2reg_lo = c2reg_lo/rU**2      # K2reg = c2reg/rho^2, gaussian e^{-K2reg z^2}
S2z = arb(0); S3z = arb(0)
sqyJ = sqyJ0
for side in [1, -1]:
    for j in range(J):
        z1up = arb(j+1)*dz; z0lo = arb(j)*dz
        if side == 1: Eb = iv(1, 1 + sqyJ*z1up)
        else:         Eb = iv(1 - sqyJ*z1up, 1)
        gw = (-K2reg_lo*z0lo**2).exp()
        Ph4u = up(Ph4reg(Eb)); G2u = up(G2reg(Eb))
        S2z += gw*(z1up**4*Ph4u/(24*lo(R)**4) + z1up**2*G2u/(2*lo(R)**2))*dz
        Ph3u = up(Ph3reg(Eb)); G1u = up(G1reg(Eb))
        Vu = Ph3u*z1up**3/(6*lo(R)**3) + G1u*z1up/lo(R)
        qe = Lip_eta*sqyJ*z1up**3/lo(R)**2
        if side == -1: qe = qe + G1u*sqyJ*z1up/lo(R)
        S3z += gw*(Vu**2/2)*up(qe.exp())*dz
# tails z in [zstar, zmax] subset [zstar, oo)
Ph4G = up(Ph4reg(ERA)); G2G = up(G2reg(ERA)); Ph3G = up(Ph3reg(ERA)); G1G = up(G1reg(ERA))
K = K2reg_lo; eK = (-K*zstar**2).exp()
I2 = eK*(zstar/(2*K) + 1/(4*K**2*zstar))
I4 = eK*(zstar**3/(2*K) + 3*zstar/(4*K**2) + 3/(8*K**3*zstar))
I6 = eK*(zstar**5/(2*K) + 5*zstar**3/(4*K**2) + 15*zstar/(8*K**3) + 15/(16*K**4*zstar))
S2z += Ph4G/(24*lo(R)**4)*I4 + G2G/(2*lo(R)**2)*I2
Aq = Ph3G/(6*lo(R)**3); Bq = G1G/lo(R)
# q over the tail z in [zstar, zmax]:  sqrt(y) z^3 <= sqrt(y_J) w_J^3 rho_hat^3/(2 c2reg_lo)^{3/2}
#   (since z <= zmax = w rho/sqrt(2 c2reg) and sqrt(y) w^3 <= sqrt(y_J) w_J^3 by monotonicity)
q_tail_exp = Lip_eta*sqyJ*w_J**3*rU**3/(2*c2reg_lo)**arb("1.5")/lo(R)**2 \
           + G1G*sqyJ*w_J*rU/(2*c2reg_lo).sqrt()/lo(R)
q_up = q_tail_exp.exp()
S3z += up(q_up)/2*(Aq**2*I6 + 2*Aq*Bq*I4 + Bq**2*I2)
W2reg = rU**arb("1.5")*P0reg*S2z
W3reg = rU**arb("1.5")*P0reg*S3z
print(f"S2z={S2z}, S3z={S3z}")
print(f"W2reg = {W2reg} -> C-part {float((pi/2)**arb('1.5')*arb(W2reg)):.4f}")
print(f"W3reg = {W3reg} -> C-part {float((pi/2)**arb('1.5')*arb(W3reg)):.4f}")

# ---------- TM ----------
p2reg = (pi/2)*R**3*(D2_f(T).sqrt())/(T**2+rt2*T+1)   # p2mag = y*p2reg
p2reg_lo, p2reg_up = lo(p2reg), up(p2reg)
denreg_lo = p2reg_lo*w_J/(2*c2reg_up).sqrt() - omreg_up*sqyJ
assert denreg_lo > 0, "TM denreg lost sign"
ecw = (rt2/pi).sqrt()*(-c_w/2).exp()   # sqrt(2/pi) e^{-c_w/2}
TMreg = 2*rU**arb("2.5")*P0reg*ecw*(1/denreg_lo +
        p2reg_up*(2*c2reg_up).sqrt()/(denreg_lo**2*2*c2reg_lo*w_J))
print(f"TMreg = {TMreg} -> C-part {float((pi/2)**arb('1.5')*arb(TMreg)):.4f}")

# ---------- TR, TLm ----------
GRreg = sqpi/4*2**arb("2.25")          # GR  <= t^{5/2} GRreg   (2^{9/4})
GLm2reg = sqpi/4*2**arb("4.75")        # GLm2<= t^{5/2} GLm2reg (2^{19/4})
TRreg = rU**arb("2.5")*GRreg*eRe0_up*ecw*up(E1R)/(w_J*lamRreg)
TLmreg = rU**arb("2.5")*GLm2reg*eRe0_up*ecw*up(E1L)/(w_J*lamLreg)
print(f"TRreg = {TRreg} -> C-part {float((pi/2)**arb('1.5')*arb(TRreg)):.4f}")
print(f"TLmreg = {TLmreg} -> C-part {float((pi/2)**arb('1.5')*arb(TLmreg)):.4f}")

# ---------- TLf (direct N^{3/2} bound) ----------
eta2 = iv(arb("0.5"), arb("0.5"))
lam2 = -a + (pi*lo(R)/2)*lo(Gh(eta2, T, R))
assert lam2 > 0, "lam2 lost sign"
M1 = 1 - eta2/(2*R)
ureg = rt2*R/eta2 + R**2*Y/eta2**2
u = Y*ureg
logcf = iv(lo((1+up(u)).log()/up(u)), 1)     # log(1+u)/u decreasing, =1 at u=0
L2 = (rt2/4)*ureg*logcf
MmL = M1 - L2
assert up(MmL) < 0, "M1-L2 not < 0"
K_J = pi*(-up(MmL))/rt2
assert y_J < K_J/3, "y^{-3}e^{-K/y} not increasing on box"
G1h = up(sqpi/4*(1/rt2/2).gamma()/(1/rt2/2+arb("1.5")).gamma())
TLfN = G1h*(-K_J/y_J).exp()/(y_J**3*lam2)
print(f"MmL_up = {up(MmL)}, lam2 = {lam2}, TLf*N^1.5 = {TLfN}")

# ---------- TL0 (direct N^{3/2} bound) ----------
cc0 = (2+rt2)/4
c0 = (1/(2*cc0.sqrt())).log()
ex0J = c0 + (pi/N_J).sqrt()
LpN = c0 + (pi/N_J).sqrt()/2 + arb("1.5")/N_J - 1/(N_J-1)
assert LpN < 0, "TL0 N^{3/2} not decreasing at N_J"
TL0N = rt2*(ex0J*N_J).exp()*N_J**arb("1.5")/(N_J-1)
print(f"TL0*N^1.5 = {TL0N}")

# ---------- D factor ----------
Ebd_up = Tval*(arb("0.75")+rt2/3) + (arb(459)/24)*Tval**3/(1-3*Tval)
d1 = y_J
d2 = (rU*(arb("0.75")+rt2/3) + (arb(459)/24)*Tval**2*rU/(1-3*Tval))*Ebd_up.exp()
d3 = rU*(6 + 5*Tval + 2*Tval**2)/4
d4 = mSreg_up*(y_J*mSreg_up).exp()
invp2reg = (2/(pi*R**3))*((numR_f(T)**2 + (T+rt2)**2).sqrt())/D2_f(T)
d5 = omreg_up**2*up(invp2reg)/2*(y_J*omreg_up**2*up(invp2reg)/2).exp()  # |e^x-1|<=|x|e^|x|, x=om^2 invp2/2
zz = rt2*R*Y/(2+rt2*R*Y); z2 = zz*zz
A1mreg = iv(0, up(2*R**2*Y/((2+rt2*R*Y)**2*(1+z2))))
breg = -SG/(2*R) + (-2*Y*SG + 2*rt2*SG + rt2*R - 2*rt2*R*A1mreg)/(2*rt2*(2+rt2*R*Y))
breg_abs = max(abs(lo(breg)), abs(up(breg)))
d6 = (pi/rt2)*breg_abs
Sd = d1+d2+d3+d4+d5+d6
Dreg = (pi/rt2)*Sd*(y_J*Sd).exp()
print(f"d1..d6 = {d1}, {d2}, {d3}, {d4}, {d5}, {d6}")
print(f"Dreg = {Dreg}")

# ---------- assembly ----------
EIreg = arb(W1reg)+arb(W2reg)+arb(W3reg)+arb(TMreg)+arb(TRreg)+arb(TLmreg)
C = (2/pi)*((pi/2)**arb("1.5")*EIreg + arb(TLfN) + arb(TL0N)) + arb(Dreg)
print()
print(f"=== FINAL CELL y in [0, {float(y_J):.6e}] (N >= 1e10): C = {C.upper()} ===")
print(f"    ({time.time()-t0:.1f}s)")
import json
json.dump({"C": str(C.upper()), "W1": str(W1reg), "W2": str(W2reg), "W3": str(W3reg),
           "TM": str(TMreg), "TR": str(TRreg), "TLm": str(TLmreg),
           "TLfN": str(TLfN), "TL0N": str(TL0N), "D": str(Dreg)}, open("final_cell_J192.json","w"), indent=1)
with open('progress.log','a') as f:
    f.write(f"turn3: FINAL CELL certified: C = {float(C.upper()):.4f} for all N>=1e10 (arb). Pieces printed.\n")
