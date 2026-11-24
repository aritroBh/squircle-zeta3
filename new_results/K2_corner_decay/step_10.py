
import mpmath as mp
mp.mp.dps = 60

def gamma_coeffs(c, M1):
    Q = [mp.mpf(1), -2*c, mp.mpf(1)] + [mp.mpf(0)]*(M1-2)
    b = [mp.sqrt(Q[0])]
    for n in range(1, M1+1):
        s = Q[n] - mp.fsum(b[k]*b[n-k] for k in range(1, n))
        b.append(s/(2*b[0]))
    return b

alpha = mp.pi/2
c = mp.cos(alpha)
gam = gamma_coeffs(c, 12)
M1, J1 = 6, 6

def S_func(g, u):
    beta = [mp.mpc(1)]
    for j in range(J1):
        beta.append(beta[-1]*(-(1+u+j))/(j+1))
    A = mp.exp(mp.log(g)*(u-1)); B = mp.exp(-mp.log(g)*(u-1))
    C = mp.exp(-mp.log1p(g)*(1+u))
    Q = 1 - 2*c*B + B*B
    F = (C*mp.sqrt(Q) - 1)*A
    S = F
    Bp = mp.mpc(1)
    for m in range(1, M1+1):
        S -= gam[m]*Bp; Bp *= B
    Ap = A*g
    for j in range(1, J1+1):
        S -= beta[j]*Ap; Ap *= g
    return S

u = mp.mpc(0)
ana = mp.mpc(0)
for m in range(1, M1+1): ana += gam[m]/(m-(m-1)*u)
beta = [mp.mpc(1)]
for j in range(J1): beta.append(beta[-1]*(-(1+u+j))/(j+1))
for j in range(1, J1+1): ana += beta[j]/(u+j)

v_quad = mp.quad(lambda g: S_func(g,u), [0,1]) + ana
v_direct = mp.quad(lambda g: ((1+g)**(-1-u)*mp.sqrt(1+g**(2-2*u)-2*c*g**(1-u))-1)*g**(u-1), [0,1])
print("quad(S)+ana =", mp.nstr(v_quad, 25))
print("direct      =", mp.nstr(v_direct, 25))
print("diff        =", mp.nstr(abs(v_quad-v_direct), 5))

for h,J in [(mp.mpf('0.16'),34),(mp.mpf('0.08'),68),(mp.mpf('0.04'),136)]:
    tot = mp.mpc(0)
    for jj in range(-J, J+1):
        t = jj*h; st = mp.sinh(t)
        g = mp.tanh(mp.pi/2*st)
        w = (mp.pi/2)*mp.cosh(t)/(mp.cosh(mp.pi/2*st)**2)
        if g<=0 or g>=1: continue
        tot += w*S_func(g,u)
    print("h=",mp.nstr(h,3),"TS(S)+ana =", mp.nstr(tot*h+ana, 25), " err=", mp.nstr(abs(tot*h+ana-v_direct),3))
