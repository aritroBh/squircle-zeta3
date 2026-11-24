from flint import arb, acb, ctx
ctx.dps = 70
pi = arb.pi(); rt2 = arb(2).sqrt(); sqpi = pi.sqrt()
a = pi/(2*rt2)
def iv(lo, hi):
    lo = arb(lo).lower(); hi = arb(hi).upper()
    return arb((lo+hi)/2, (hi-lo)/2)
def up(x): return arb(x.upper())
def lo(x): return arb(x.lower())
N_J = arb(10)**10
y_J = (pi/(2*N_J)).sqrt()
Y  = iv(0, y_J); SG = iv(0, 1); R  = iv(1, 1 + y_J**2)
rU = up(R); Tval = up(R*y_J); T  = iv(0, Tval)
c_w = arb("1.0")
L_J = (-y_J.log()); al_ = (pi/2).log() + c_w
w_J = (al_ + 2*L_J).sqrt()
sqyJ0 = y_J.sqrt()
def numR_f(Tb): return Tb**3 + 2*rt2*Tb**2 + 4*Tb + rt2
def D2_f(Tb):  return Tb**2 + 2*rt2*Tb + 4
c2reg = (pi*R**3/4)*numR_f(T)/(1+rt2*T+T**2)**2
c2reg_lo, c2reg_up = lo(c2reg), up(c2reg)
Omn = rt2*R + R**2*Y - 2*Y*SG - Y**3*SG**2
omreg = (pi/(2*rt2))*Omn/(1+rt2*R*Y+R**2*Y**2)
omreg_up = up(omreg)
mSreg_up = up(pi*R**3/4)*rt2/3
eRe0_up = (y_J*mSreg_up).exp()
G0reg = sqpi/4 * 2**arb("1.5")
P0reg = G0reg * eRe0_up
p2reg = (pi/2)*R**3*(D2_f(T).sqrt())/(T**2+rt2*T+1)
p2reg_lo, p2reg_up = lo(p2reg), up(p2reg)
denreg_lo = p2reg_lo*w_J/(2*c2reg_up).sqrt() - omreg_up*sqyJ0
ecw = (rt2/pi).sqrt()*(-c_w/2).exp()
TMreg = 2*rU**arb("2.5")*P0reg*ecw*(1/denreg_lo + p2reg_up*(2*c2reg_up).sqrt()/(denreg_lo**2*2*c2reg_lo*w_J))
print("denreg_lo =", denreg_lo)
print("2*rU^2.5*P0reg*ecw =", 2*rU**arb("2.5")*P0reg*ecw)
print("1/denreg_lo =", 1/denreg_lo)
print("second =", p2reg_up*(2*c2reg_up).sqrt()/(denreg_lo**2*2*c2reg_lo*w_J))
print("TMreg =", TMreg)
print("first part of TMreg =", 2*rU**arb("2.5")*P0reg*ecw/denreg_lo)
# compare with my mpmath recomputation: denreg_lo ≈ 10.32278, 2*1*1.25332*0.483941 = 1.212933
# => TMreg_first should be 1.212933/10.32278 = 0.117502 ; code prints?
