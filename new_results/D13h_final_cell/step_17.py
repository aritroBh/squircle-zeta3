# step_14.py - validate W3: my z-sum (sharp, at point) vs step_38's segmented sum vs true integral
from flint import arb, acb, ctx
ctx.dps = 60
pi = arb.pi(); rt2 = arb(2).sqrt(); sqpi = pi.sqrt()
def Hc(s): return sqpi/4*(s/2).gamma()/(s/2+arb("1.5")).gamma()
b = (acb(0,1)*pi/4).exp()
def gF(x): return Hc(b*x)/x
def Phi3b(x,N): return 2*N*(1+3*x+3*x**2)/(x**3*(1+rt2*x+x*x)**arb("1.5"))
def G1bound(x):
    t = x/(2*rt2); return arb("0.75")*(1/t**2+1/t)+1/x
def curv2(y, rho, eta):
    num = rt2*eta**3*rho + 4*eta**2*rho**2*y + 2*rt2*eta*rho**3*y**2 + rho**4*y**3
    den = eta**2*(eta**2+rt2*eta*rho*y+rho**2*y**2)**2
    return (pi*rho**2*y/4)*num/den
def h3(T): return rt2/3/(1+rt2*T+T*T)

N = arb(10)**9
y = (pi/(2*N)).sqrt(); sig = arb("0.3"); rho = 1+y*y*sig; t = rho*y; x0 = 1/t
numR = t**3+2*rt2*t**2+4*t+rt2
c2 = (pi/4)*rho**3*y*numR/(1+rt2*t+t*t)**2
c_w = arb(1); w = (N.log()+c_w).sqrt(); d = w/(2*c2).sqrt()
mS = (N*t**3/2)*h3(t)
P0t = abs(gF(x0))*mS.exp()

# (i) TRUE integral: q with sharp curvature variation along path
M = 3000
W3true = arb(0)
for side in [1,-1]:
    du = d/M
    for j in range(M):
        u = (j+arb("0.5"))*du
        x = x0 + side*u
        eta = x/x0
        Rb = Phi3b(x,N)*u**3/6
        Gb = G1bound(x)*u
        V = Rb+Gb
        c2eta = curv2(y, rho, eta)
        qe = (c2-c2eta)*u*u
        if side == -1: qe = qe + G1bound(x)*u
        W3true += (-c2*u*u).exp()*(V*V/2)*qe.exp()*du
W3true *= P0t
print("W3 true*N^1.5      =", float(W3true*N**arb("1.5")))

# (ii) step_38 segmented (J=24) at point
J = 24
W338 = arb(0)
for side in [1,-1]:
    for j in range(J):
        u0 = d*j/J; u1 = d*(j+1)/J
        xsegL = x0 if side==1 else x0-u1
        M3 = Phi3b(xsegL,N); G1 = G1bound(xsegL)
        Rb = M3*u1**3/6; Gb = G1*u1; V = Rb+Gb
        c2side = curv2(y, rho, (x0+side*d)/x0 if side==1 else (x0-d)/x0)
        qe = (c2-c2side)*u1**2
        if side == -1: qe = qe + G1*u1
        gw = (-c2*u0**2).exp()
        W338 += P0t*gw*(V*V/2)*max(arb(1), qe.exp())*(d/J)
print("W3 step38*N^1.5    =", float(W338*N**arb("1.5")))

# (iii) my z-space sum (J=48, zstar=3.0) at point + tail
c2reg = c2/y; K2reg = c2reg/rho**2
P0reg = P0t/t**arb("2.5")
def Ph3regP(eta):
    return pi*rho**4*(rho**2*y**2+3*eta*rho*y+3*eta**2)/(eta**3*(eta**2+rt2*eta*rho*y+rho**2*y**2)**arb("1.5"))
def G1regP(eta):
    return 6*rho**2*y/eta**2 + (arb("1.5")*rt2+1)*rho/eta
sqy = y**arb("0.5")
# Lip_eta pointwise: |d/deta curv2reg| at eta=1 (dominant), computed numerically
eps = arb("1e-6")
def c2regP(eta):
    return curv2(y, rho, eta)/y
Lip_pt = abs((c2regP(1+eps)-c2regP(1-eps))/(2*eps))
J = 48; zstar = arb("3.0"); dz = zstar/J
S3z = arb(0)
for side in [1,-1]:
    for j in range(J):
        z1 = (j+1)*dz; z0 = j*dz
        eta = 1 + side*z1*sqy
        gw = (-K2reg*z0**2).exp()
        Vu = Ph3regP(eta)*z1**3/(6*rho**3) + G1regP(eta)*z1/rho
        qe = Lip_pt*sqy*z1**3/rho**2
        if side == -1: qe = qe + G1regP(eta)*sqy*z1/rho
        S3z += gw*(Vu**2/2)*qe.exp()*dz
# tail
K = K2reg; eK = (-K*zstar**2).exp()
I2 = eK*(zstar/(2*K)+1/(4*K**2*zstar)); I4 = eK*(zstar**3/(2*K)+3*zstar/(4*K**2)+3/(8*K**3*zstar)); I6 = eK*(zstar**5/(2*K)+5*zstar**3/(4*K**2)+15*zstar/(8*K**3)+15/(16*K**4*zstar))
Aq = Ph3regP(1)/(6*rho**3); Bq = G1regP(1)/rho
q_up = (Lip_pt*sqy*(2*zstar)**3/rho**2 + G1regP(1)*sqy*(2*zstar)/rho).exp()
S3z += q_up/2*(Aq**2*I6+2*Aq*Bq*I4+Bq**2*I2)
W3mine = (pi/2)**arb("1.5")*rho**arb("1.5")*P0reg*S3z
print("W3 mysum (pi/2)^1.5*rho^1.5 P0reg S3z =", float(W3mine))
print("W3 mine/y^3*(2/pi) -> C part:", float((2/pi)*W3mine))
print("my final-cell W3 C-part was 43.77; step38 point gives 47.6; true =", float(W3true*N**arb("1.5")))
