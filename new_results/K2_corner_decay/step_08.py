
# step_02.py --- prototype: M(u;alpha) at complex u via subtracted tanh-sinh; verify
import mpmath as mp

mp.mp.dps = 60

def gamma_coeffs(c, M1):
    # sqrt(1 - 2 c z + z^2) = sum gamma_m z^m
    Q = [mp.mpf(1), -2*c, mp.mpf(1)] + [mp.mpf(0)]*(M1-2)
    b = [mp.sqrt(Q[0])]
    for n in range(1, M1+1):
        s = Q[n] - mp.fsum(b[k]*b[n-k] for k in range(1, n))
        b.append(s/(2*b[0]))
    return b  # gamma_0..gamma_M1

def make_M_func(alpha, M1=6, J1=6):
    c = mp.cos(alpha)
    gam = gamma_coeffs(c, M1)   # gamma_0..gamma_M1
    s2a2 = mp.sin(alpha)**2/2
    def Mval(u):
        # subtraction terms: m=1..M1 of gam_m B^{m-1}, j=1..J1 of beta_j A g^j
        beta = [mp.mpc(1)]
        for j in range(J1):
            beta.append(beta[-1]*(-(1+u+j))/(j+1))
        # analytic part
        ana = mp.mpc(0)
        for m in range(1, M1+1):
            ana += gam[m]/(m-(m-1)*u)
        for j in range(1, J1+1):
            ana += beta[j]/(u+j)
        # tanh-sinh on g in (0,1)
        h = mp.mpf('0.16'); J = 34
        tot = mp.mpc(0)
        for jj in range(-J, J+1):
            t = jj*h
            st = mp.sinh(t)
            g = mp.tanh(mp.pi/2*st)
            w = (mp.pi/2)*mp.cosh(t)/(mp.cosh(mp.pi/2*st)**2)
            if g <= 0 or g >= 1:
                continue
            A = mp.exp(mp.log(g)*(u-1))     # g^{u-1}
            B = mp.exp(-mp.log(g)*(u-1))    # g^{1-u}
            C = mp.exp(-mp.log1p(g)*(1+u))  # (1+g)^{-1-u}
            Q = 1 - 2*c*B + B*B
            F = (C*mp.sqrt(Q) - 1)*A
            S = F
            Bp = mp.mpc(1)
            for m in range(1, M1+1):
                S -= gam[m]*Bp
                Bp *= B
            Ap = A*g
            for j in range(1, J1+1):
                S -= beta[j]*Ap
                Ap *= g
            tot += w*S
        return tot*h + ana
    return Mval

# verification 1: M(0; alpha) vs direct integral
for alpha in [mp.pi/2, mp.pi/3, mp.mpf('2')*mp.pi/3, mp.mpf('3')*mp.pi/4, mp.mpf('5')*mp.pi/6]:
    Mf = make_M_func(alpha)
    v1 = Mf(mp.mpc(0))
    c = mp.cos(alpha)
    v2 = mp.quad(lambda g: ((1+g)**(-1)*mp.sqrt(1+g*g-2*c*g)-1)/g, [0,1])
    print("alpha=%s  M(0) sub=%s  direct=%s  diff=%.2e" % (
        mp.nstr(alpha/mp.pi,5), mp.nstr(v1.real,20), mp.nstr(v2,20), abs(v1-v2)))
