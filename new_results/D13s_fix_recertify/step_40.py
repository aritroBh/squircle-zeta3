
import os
WORK = 'research/08_theory/D13s_fix_recertify/D13s_work'
FIX = os.path.join(WORK, 'fixed')
src = open(os.path.join(WORK, 'D13h_final_cell__step_10_J192.py')).read()

# ---- B1 fix ----
old = """ecw = (rt2/pi).sqrt()*(-c_w/2).exp()   # sqrt(2/pi) e^{-c_w/2}"""
new = """ecw = (arb(2)/pi).sqrt()*(-c_w/2).exp()   # sqrt(2/pi) e^{-c_w/2}  [D13s B1 FIX: was (rt2/pi).sqrt() = sqrt(sqrt(2)/pi), off by 2^{-1/4}]"""
assert old in src; src = src.replace(old, new)

# ---- H2 fix: proved TL0 ----
old = """cc0 = (2+rt2)/4
c0 = (1/(2*cc0.sqrt())).log()
ex0J = c0 + (pi/N_J).sqrt()
LpN = c0 + (pi/N_J).sqrt()/2 + arb("1.5")/N_J - 1/(N_J-1)
assert LpN < 0, "TL0 N^{3/2} not decreasing at N_J"
TL0N = rt2*(ex0J*N_J).exp()*N_J**arb("1.5")/(N_J-1)"""
new = """# D13s H2 FIX (PROVED bound, see fixed step_38.py comment):
# I0 <= (2+rt2)^{-N/2} * 2/(N-2); TL0*N^{3/2} <= 2 (2+rt2)^{-N/2} e^{theta} N^{3/2}/(N-2), theta=sqrt(pi N).
# decreasing in N for N >= N_J: d/dN log[...] = -c0 + sqrt(pi)/(2 sqrt N) + 1.5/N - 1/(N-2) < 0,
# and that log-derivative is itself decreasing in N (so stays < 0 for all N >= N_J).
c0 = (2+rt2).log()/2
LpN = -c0 + (pi/N_J).sqrt()/2 + arb("1.5")/N_J - 1/(N_J-2)
assert LpN < 0, "TL0 N^{3/2} not decreasing at N_J"
TL0N = arb(2)*(-c0*N_J + (pi*N_J).sqrt()).exp()*N_J**arb("1.5")/(N_J-2)"""
assert old in src; src = src.replace(old, new)

# ---- TLfN (pi/2)^{3/2} factor fix ----
old = """TLfN = G1h*(-K_J/y_J).exp()/(y_J**3*lam2)"""
new = """TLfN = G1h*(pi/2)**arb("1.5")*(-K_J/y_J).exp()/(y_J**3*lam2)   # [D13s: (pi/2)^{3/2} factor restored (review's 'smaller item'); value ~1e-19229, harmless]"""
assert old in src; src = src.replace(old, new)

# ---- H1 fix: exact Btrue regular form replaces breg ----
old = """zz = rt2*R*Y/(2+rt2*R*Y); z2 = zz*zz
A1mreg = iv(0, up(2*R**2*Y/((2+rt2*R*Y)**2*(1+z2))))
breg = -SG/(2*R) + (-2*Y*SG + 2*rt2*SG + rt2*R - 2*rt2*R*A1mreg)/(2*rt2*(2+rt2*R*Y))
breg_abs = max(abs(lo(breg)), abs(up(breg)))
d6 = (pi/rt2)*breg_abs"""
new = """zz = rt2*R*Y/(2+rt2*R*Y); z2 = zz*zz
A1mreg = iv(0, up(2*R**2*Y/((2+rt2*R*Y)**2*(1+z2))))
breg = -SG/(2*R) + (-2*Y*SG + 2*rt2*SG + rt2*R - 2*rt2*R*A1mreg)/(2*rt2*(2+rt2*R*Y))
breg_abs = max(abs(lo(breg)), abs(up(breg)))
# D13s H1 FIX: replace breg by the PROVED EXACT regular form of Btrue = (phi+pi/4)/((pi/rt2) y)
# on the saddle branch rho(t) = sqrt(1+t^2/(1+rt2 t)), y = t/rho (saddle eq RePhi'(1/t)=0 PROVED):
#   phi+pi/4 = a(1-2rho)/t + N atan(zz) + pi/4  [arg(1+b x0) = atan(1/(1+rt2 t)); pi/4 - that = atan(zz)]
#   Btrue = rho*B1/(2+rt2 t),  B1 = 1/2 + rt2 s1 + s1^2 - rho^2 A1n,
#   s1 = (1-rho)/t = -t/((1+rt2 t)(1+rho))   [rationalized, regular],
#   A1n = A1'/t^2, A1' = 1 - atan(zz)/zz in [zz2/3 - zz2^2/5, zz2/3] (alternating series, |zz|<=1),
#   so A1n in [2/(3(2+rt2 t)^2) - 4 t^2/(5(2+rt2 t)^4), 2/(3(2+rt2 t)^2)]   (regular at t=0).
# Verified vs direct Btrue: worst rel resid 8.5e-37 on t in [1e-12,1] (verify_Btrue.py, 60 dps).
Tbox = iv(0, Tval)
rhoB = (1 + Tbox**2/(1+rt2*Tbox)).sqrt()
s1B = -Tbox/((1+rt2*Tbox)*(1+rhoB))
A1nB = iv(lo(2/(3*(2+rt2*Tbox)**2) - 4*Tbox**2/(5*(2+rt2*Tbox)**4)),
          up(2/(3*(2+rt2*Tbox)**2)))
B1box = arb("0.5") + rt2*s1B + s1B**2 - rhoB**2*A1nB
Btrue_box = rhoB*B1box/(2+rt2*Tbox)
B_abs = max(abs(lo(Btrue_box)), abs(up(Btrue_box)))
print(f"Btrue box over final cell = {Btrue_box}, B_abs = {B_abs} (limit 1/6 = {arb(1)/6})")
assert breg_abs >= B_abs, "old breg_abs does not dominate B_abs"   # old breg was harmless (|Btrue|<=B_abs<=breg_abs)
d6 = (pi/rt2)*B_abs"""
assert old in src; src = src.replace(old, new)

open(os.path.join(FIX, 'step_10_J192.py'), 'w').write(src)
print("fixed/step_10_J192.py written,", len(src), "bytes")
