from flint import arb, ctx
ctx.dps = 70
pi = arb.pi(); rt2 = arb(2).sqrt()
a = pi/(2*rt2)
def iv(lo, hi):
    lo = arb(lo).lower(); hi = arb(hi).upper()
    return arb((lo+hi)/2, (hi-lo)/2)
def up(x): return arb(x.upper())
def lo(x): return arb(x.lower())

N_J = arb(10)**10
y_J = (pi/(2*N_J)).sqrt()
Y  = iv(0, y_J); SG = iv(0,1); R = iv(1, 1+y_J**2); rU = up(R)
Tval = up(R*y_J); T = iv(0, Tval)
c_w = arb("1.0"); L_J = -y_J.log(); al_ = (pi/2).log()+c_w; w_J = (al_+2*L_J).sqrt()
def numR_f(Tb): return Tb**3 + 2*rt2*Tb**2 + 4*Tb + rt2
c2reg = (pi*R**3/4)*numR_f(T)/(1+rt2*T+T**2)**2
c2reg_up = up(c2reg)
g1R = arb("1.319008055487752087349884188819088309400330572142180472831699842794441")
Q1 = R**2 - SG*(2+Y**2*SG)*(1+rt2*T)
Q1abs = max(abs(lo(Q1)), abs(up(Q1)))
y15 = up(Y**arb("1.5"))
print("Q1 =", Q1, "Q1abs =", Q1abs)
print("y_J^1.5 =", y15)
term1 = (pi*lo(R)**2/2)*g1R/(2*c2reg_up).sqrt()
term2 = a*y15*Q1abs/w_J
print("term1 =", term1)
print("term2 =", term2)
print("lamRreg =", term1 - term2)
print("w_J =", w_J, "c2reg_up =", c2reg_up)
