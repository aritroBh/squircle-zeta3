
import mpmath as mp
mp.mp.dps = 70

# --- tiny jet algebra in p, order 3 ---
def jmul(A, B):
    a0,a1,a2,a3 = A; b0,b1,b2,b3 = B
    return (a0*b0, a0*b1+a1*b0, a0*b2+2*a1*b1+a2*b0, a0*b3+3*a1*b2+3*a2*b1+a3*b0)
def jadd(A,B): return tuple(x+y for x,y in zip(A,B))
def jscale(A,c): return tuple(c*x for x in A)
def jdiv(A,B):
    a0,a1,a2,a3=A; b0,b1,b2,b3=B
    c0 = a0/b0
    c1 = (a1 - c0*b1)/b0
    c2 = (a2 - 2*c1*b1 - c0*b2)/b0
    c3 = (a3 - 3*c2*b1 - 3*c1*b2 - c0*b3)/b0
    return (c0,c1,c2,c3)
def jlog(A):
    a0,a1,a2,a3=A
    b0 = mp.log(a0); b1 = a1/a0; b2 = a2/a0 - (a1/a0)**2
    b3 = a3/a0 - 3*a1*a2/a0**2 + 2*(a1/a0)**3
    return (b0,b1,b2,b3)
def jexp(A):
    a0,a1,a2,a3=A
    b0 = mp.exp(a0); b1 = b0*a1; b2 = b0*(a2+a1**2); b3 = b0*(a3+3*a1*a2+a1**3)
    return (b0,b1,b2,b3)
def jsqrt(A):
    return jexp(jscale(jlog(A), mp.mpf('0.5')))

def term_jets(phi, m):
    # returns for each i: a_i, a_i', a_i'' (phi derivatives), with sign conventions
    out = []
    for i in range(1, m+1):
        th = phi - i*mp.pi/m
        c = mp.cos(th); s = mp.sin(th)
        a = abs(c); ap = -mp.sign(c)*s; app = -a
        out.append((a, ap, app))
    return out

def S_jet(phi, m, p):
    J = (mp.mpf(0),)*4
    for a, ap, app in term_jets(phi, m):
        la = mp.log(a) if a > 0 else mp.mpf('-1e999')
        ap_p = a**p
        J = jadd(J, (ap_p, ap_p*la, ap_p*la**2, ap_p*la**3))
    return J

def Sphi_jet(phi, m, p):
    J = (mp.mpf(0),)*4
    for a, ap, app in term_jets(phi, m):
        la = mp.log(a) if a > 0 else mp.mpf('-1e999')
        # p * a^(p-1) * ap : jet
        base = a**(p-1)*ap
        J = jadd(J, (base*p, base*(1+p*la), base*la*(2+p*la), base*la**2*(3+p*la)))
    return J

def Sphiphi_jet(phi, m, p):
    J = (mp.mpf(0),)*4
    for a, ap, app in term_jets(phi, m):
        la = mp.log(a) if a > 0 else mp.mpf('-1e999')
        # d/dphi [p a^(p-1) ap] = p(p-1)a^(p-2) ap^2 + p a^(p-1) app
        t1 = p*(p-1)*a**(p-2)*ap**2
        t2 = p*a**(p-1)*app
        base = t1 + t2
        # jet in p of base: base is built from terms; do termwise
        # t1 = (p^2-p) a^(p-2) ap^2
        J1 = (ap**2 * a**(p-2) * (p**2-p),
              ap**2 * a**(p-2) * (2*p-1 + (p**2-p)*la),
              ap**2 * a**(p-2) * (2 + 2*(2*p-1)*la + (p**2-p)*la**2),
              ap**2 * a**(p-2) * (6*la + 3*(2*p-1)*la**2 + (p**2-p)*la**3))
        # t2 = p a^(p-1) app
        J2 = (app * a**(p-1) * p,
              app * a**(p-1) * (1 + p*la),
              app * a**(p-1) * (2*la + p*la**2),
              app * a**(p-1) * (3*la**2 + p*la**3))
        J = jadd(J, jadd(J1, J2))
    return J

def speed_jet(phi, m, p):
    SJ = S_jet(phi, m, p)
    SF = Sphi_jet(phi, m, p)
    SFF = Sphiphi_jet(phi, m, p)
    pinv = jdiv((mp.mpf(1),0,0,0),(p,1,0,0))
    Lrho = jscale(jlog(SJ), None) if False else jmul((-pinv[0],-pinv[1],-pinv[2],-pinv[3]), jlog(SJ))
    rhoJ = jexp(Lrho)
    # Lphi = -(1/p) SF/S
    LphiJ = jmul((-pinv[0],-pinv[1],-pinv[2],-pinv[3]), jdiv(SF, SJ))
    # rhophi = rho * Lphi
    rhophiJ = jmul(rhoJ, LphiJ)
    # speed^2 = rho^2 + rhophi^2
    sp2 = jadd(jmul(rhoJ,rhoJ), jmul(rhophiJ,rhophiJ))
    return jsqrt(sp2)

def P_deriv(m, k):
    pts = sorted(set([(mp.pi/2 + mp.pi*i/m) % (2*mp.pi) for i in range(2*m)] +
                     [(mp.pi*i/m) % (2*mp.pi) for i in range(2*m)]))
    pts = [x for x in pts if mp.mpf('1e-12') < x < 2*mp.pi-mp.mpf('1e-12')]
    b = [mp.mpf(0)] + pts + [2*mp.pi]
    tot = mp.mpf(0)
    for a, bb in zip(b[:-1], b[1:]):
        mid = (a+bb)/2
        tot += mp.quad(lambda ph: speed_jet(ph, m, mp.mpf(1))[k], [a, mid, bb])
    return tot

print("m=2: P1,P2,P3:", [mp.nstr(P_deriv(2,k), 22) for k in [1,2,3]])
print("  target P3:", mp.nstr(6*mp.sqrt(2)*mp.zeta(3)-2*mp.sqrt(2)*mp.pi**2, 22))
print("m=3: P1,P2,P3:", [mp.nstr(P_deriv(3,k), 22) for k in [1,2,3]])
print("  paper P1 = 2 sqrt3 ln2 =", mp.nstr(2*mp.sqrt(3)*mp.log(2), 22))
print("  C3*zeta3 =", mp.nstr(-9*mp.sqrt(3)/8*mp.zeta(3), 22))
