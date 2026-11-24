
import mpmath as mp
mp.mp.dps = 60

def gamma_coeffs(c, M1):
    Q = [mp.mpf(1), -2*c, mp.mpf(1)] + [mp.mpf(0)]*(M1-2)
    b = [mp.sqrt(Q[0])]
    for n in range(1, M1+1):
        s = Q[n] - mp.fsum(b[k]*b[n-k] for k in range(1, n))
        b.append(s/(2*b[0]))
    return b

def make_M_func(alpha, M1=6, J1=6):
    c = mp.cos(alpha)
    gam = gamma_coeffs(c, M1)
    def Mval(u, h=mp.mpf('0.12'), J=40):
        beta = [mp.mpc(1)]
        for j in range(J1):
            beta.append(beta[-1]*(-(1+u+j))/(j+1))
        ana = mp.mpc(0)
        for m in range(1, M1+1): ana += gam[m]/(m-(m-1)*u)
        for j in range(1, J1+1): ana += beta[j]/(u+j)
        tot = mp.mpc(0)
        for jj in range(-J, J+1):
            t = jj*h; st = mp.sinh(t)
            th = mp.tanh(mp.pi/2*st)
            g = (1+th)/2
            w = (mp.pi/4)*mp.cosh(t)/(mp.cosh(mp.pi/2*st)**2)
            if g <= 0 or g >= 1: continue
            A = mp.exp(mp.log(g)*(u-1)); B = mp.exp(-mp.log(g)*(u-1))
            C = mp.exp(-mp.log1p(g)*(1+u))
            Q = 1 - 2*c*B + B*B
            S = (C*mp.sqrt(Q) - 1)*A
            Bp = mp.mpc(1)
            for m in range(1, M1+1):
                S -= gam[m]*Bp; Bp *= B
            Ap = A*g
            for j in range(1, J1+1):
                S -= beta[j]*Ap; Ap *= g
            tot += w*S
        return tot*h + ana
    return Mval

for alpha in [mp.pi/2, mp.pi/3]:
    Mf = make_M_func(alpha)
    c = mp.cos(alpha)
    v2 = mp.quad(lambda g: ((1+g)**(-1)*mp.sqrt(1+g*g-2*c*g)-1)/g, [0,1])
    for h,J in [(mp.mpf('0.12'),40),(mp.mpf('0.08'),60),(mp.mpf('0.05'),96)]:
        v1 = Mf(mp.mpc(0), h, J)
        print("alpha=%s h=%s err=%.2e" % (mp.nstr(alpha/mp.pi,4), mp.nstr(h,3), abs(v1-v2)))
    # complex u
    v2c = mp.quad(lambda g: ((1+g)**(-1-mp.mpc('0.3','0.7'))*mp.sqrt(1+g**(2-2*mp.mpc('0.3','0.7'))-2*c*g**(1-mp.mpc('0.3','0.7')))-1)*g**(mp.mpc('0.3','0.7')-1), [0,1])
    v1c = Mf(mp.mpc('0.3','0.7'), mp.mpf('0.05'), 96)
    print("   complex u=0.3+0.7i: ours=%s quad=%s err=%.2e" % (mp.nstr(v1c,15), mp.nstr(v2c,15), abs(v1c-v2c)))
