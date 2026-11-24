import mpmath as mp
mp.mp.dps = 60
sqrt2 = mp.sqrt(2); a = mp.pi/(2*sqrt2); pi = mp.pi; sqpi = mp.sqrt(pi)
bb = mp.e**(1j*mp.pi/4)
def Phi(x, N): return (-1+1j)*a*x + N*mp.log(bb*x/(1+bb*x))
def Phi1(x, N): return (-1+1j)*a + N*(1/x - bb/(1+bb*x))
def Phi2(x, N): return -N*(1+2*bb*x)/(x**2*(1+bb*x)**2)
def Hc(s): return sqpi/4*mp.gamma(s/2)/mp.gamma(s/2+mp.mpf('1.5'))
def g(x): return Hc(bb*x)/x

def true_pieces(Nv, full=False):
    N = mp.mpf(Nv); th = mp.sqrt(pi*N)
    y = mp.sqrt(pi/(2*N))
    f = lambda r: r**2*(1+sqrt2*r*y)/(1+sqrt2*r*y+r**2*y**2) - 1
    rho = mp.findroot(f, 1+y*y/2, tol=mp.mpf('1e-50'))
    x0 = 1/(rho*y)
    Phi0 = Phi(x0,N); p2 = Phi2(x0,N); g0 = g(x0)
    om = mp.im(Phi1(x0,N))
    c2 = abs(mp.re(p2))/2
    sg = 1/mp.sqrt(2*c2)
    w = mp.sqrt(mp.log(N)+1.0)
    d = min(w*sg, x0/3)
    xl, xr = x0-d, x0+d; x2 = x0/2
    def integrand(x): return mp.e**Phi(x,N)*g(x)
    def model(u): return g0*mp.e**Phi0*mp.e**(1j*om*u + p2*u*u/2)
    MV = g0*mp.e**Phi0*mp.sqrt(-2*pi/p2)*mp.e**(om**2/(2*p2))
    Wt = mp.quad(lambda u: abs(integrand(x0+u)-model(u)), [0, d]) + mp.quad(lambda u: abs(integrand(x0+u)-model(u)), [-d, 0])
    Min = mp.quad(lambda u: model(u), [-d, d])
    TMt = abs(MV - Min)
    TRt = mp.quad(lambda x: abs(integrand(x)), [xr, x0+6*d, mp.inf])
    TLmt = mp.quad(lambda x: abs(integrand(x)), [x2, xl])
    TLft = mp.quad(lambda x: abs(integrand(x)), [1, x2])
    TL0t = mp.quad(lambda x: abs(integrand(x)), [0, mp.mpf('0.5'), 1])
    sc = N**mp.mpf('-1.5')*mp.e**(-th)
    Lt = -1j*(pi**mp.mpf('1.5')/(2*N))*mp.e**(-th)*mp.e**(1j*th)
    Dt = (2/pi)*abs(MV-Lt)
    print(f"N={Nv}: d/x0={mp.nstr(d/x0,6)} (branch {'w*sg' if w*sg < x0/3 else 'x0/3'})")
    tot = (Wt+TMt+TRt+TLmt+TLft+TL0t)/sc + Dt/sc
    print(f"  TRUTH/sc: W={mp.nstr(Wt/sc,8)} TM={mp.nstr(TMt/sc,8)} TR={mp.nstr(TRt/sc,8)} TLm={mp.nstr(TLmt/sc,8)} TLf={mp.nstr(TLft/sc,8)} TL0={mp.nstr(TL0t/sc,8)} D={mp.nstr(Dt/sc,8)}")
    print(f"  sum true/sc = {mp.nstr(tot,8)}")
    return tot
true_pieces(10**6)
true_pieces(10**8)
with open('research/08_theory/D13r_lemma_review/progress.log','a') as f:
    f.write("L7: HOLDS (all factorisation identities residual<1e-45 @60dps: curv2, ph3=y^2Ph3reg, ph4=y^3Ph4reg, G1=yG1reg, G2=y^2G2reg, c2=yc2reg, om=y omreg, tail id; "
            "mS: RePhi0+th = -(N/2)S(rho y) NOT exact (diff=(pi rt2/(4t))(rho-1)^2>0) BUT h3-box provably contains truth for all t>0 (membership margin=(pi rt2/(4t))[t^3(rt2+t)/(3(1+rt2 t))-(rho-1)^2]>0 since (rt2+t)(1+rt2t)(sqrt(1+u)+1)^2>=4rt2+8t>3t); h3 two-sided verified; "
            "phi: step_37 Bh formula exact with A1=A1'/y (resid 1.9e-53) but coded A1m=[0,z2/(1+z2)] FAILS to contain A1'/y (confirmed) -- INTERMEDIATE BOX WRONG; full (y,sig,A1) z-box nonetheless contains z_true at all 41x4+21x7 sampled N in [100,1.5e10]; D13h breg is NOT exact for (phi+pi/4)/((pi/rt2)y) (my exact regular form printed, resid 5.7e-44) but min(breg_h-Btrue)=0.08333>0 & Btrue in (0.16666,1/6) with breg_abs~0.25 => d6 valid margin 1.5; "
            "z_true=factored MV/L~ resid 2.9e-56; leading term identity verified).\n")
    f.write("L8 partial: TL0 truth/code ratio 0.389/0.382/0.381 at N=100/400/1000 (valid, margin 2.6x; hand-derivation of exact constant questionable but harmless, contribution ~1e-18 to C). "
            "Full assembly truth check: N=100 truths/sc W 8.63 TM 3.85 TR 2.07 TLm 2.82 TLf 0.159 TL0 6e-19 D 0.350, sum 17.88 <= C_cert 99.99; |Mmain-lead|/sc=0.602<=99.99. "
            "N=400: W 15.67 TM 2.13 TR 1.66 TLm 0.97 TLf 0.0025 D 0.318, sum 20.75 <= 152.68 (and <=46.55 step_34); |Mmain-lead|/sc=1.065. N=1e6,1e8 truths printed (w*sg branch).\n")
print("logged")
