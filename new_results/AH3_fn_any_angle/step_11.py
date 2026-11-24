
import mpmath as mp
mp.mp.dps = 50

def phi(s, th):
    c = mp.cos(th)
    return 1 - (s + c)/mp.sqrt(1 + 2*c*s + s**2)

def S_of_U(u, th):
    C = mp.cos(th/2)**2; Z = mp.sin(th/2)**2
    return C/u - mp.cos(th) - Z*u

def legval(nu, x):
    return mp.legenp(nu, 0, x)

def Lambda(nu, th):
    # int_{cos th}^1 P_nu(t) dt = 2 sin^2(th/2) 2F1(-nu,1+nu;2;sin^2(th/2))
    z = mp.sin(th/2)**2
    return 2*z*mp.hyper([-nu, 1+nu], [2], z)

# 1) Verify E[S_theta^w] = (pi w / sin(pi w)) * Lambda_w(theta)/(1-cos theta), two ways:
print("Mellin identity checks:")
for th in [mp.pi/3, mp.pi/2, 2*mp.pi/3, mp.pi/6]:
    c = mp.cos(th)
    for w in [mp.mpf('0.3'), mp.mpf('-0.4'), mp.mpf('0.77')]:
        lhs = mp.quad(lambda s: s**(w-1)*phi(s,th), [0, mp.inf])/(1-c)
        lhs2 = mp.quad(lambda u: S_of_U(u, th)**w, [0,1])           # from explicit law
        rhs = (mp.pi*w/mp.sin(mp.pi*w))*Lambda(w, th)/(1-c)
        print(f" th={mp.nstr(th,6)} w={w}: int={mp.nstr(lhs,12)} explaw={mp.nstr(lhs2,12)} leg={mp.nstr(rhs,12)} maxdiff={mp.nstr(max(abs(lhs-rhs),abs(lhs2-rhs)),2)}")

# 2) Verify Lambda integral definition vs hypergeometric
for th in [mp.pi/3, 2*mp.pi/3]:
    for nu in [mp.mpf('0.25'), mp.mpf('0.8')]:
        i1 = mp.quad(lambda t: legval(nu, t), [mp.cos(th), 1])
        print('Lambda quad vs 2F1 diff:', mp.nstr(abs(i1 - Lambda(nu, th)), 2))
