from flint import arb, acb, ctx
import json, sys
ctx.dps = 70
pi = arb.pi(); rt2 = arb(2).sqrt(); a = pi/(2*rt2); sqpi = pi.sqrt()
b = (acb(0,1)*pi/4).exp()          # e^{i pi/4}
I1 = acb(0,1)

def iv(lo, hi):
    lo = arb(lo).lower(); hi = arb(hi).upper()
    return arb((lo+hi)/2, (hi-lo)/2)

def Hr(t):   # H on positive reals (arb)
    return sqpi/4*(t/2).gamma()/(t/2+arb("1.5")).gamma()
def Hc(s):   # H on complex (acb)
    return sqpi/4*(s/2).gamma()/(s/2+arb("1.5")).gamma()
def PhiF(x, N):                     # x arb box, N arb -> acb
    s = b*x
    return I1*pi*s/2 + N*(s/(1+s)).log()
def RePhi(x,N): return -a*x + N*x.log() - N*(1+rt2*x+x*x).log()/2
def ReP1(x,N):  return -a + N/x - N*(x+1/rt2)/(1+rt2*x+x*x)
def ReP2(x,N):  return -N/x**2 + N*x*(x+rt2)/(1+rt2*x+x*x)**2
def Phi3F(x,N):
    s = b*x
    return 2*N/x**3 - 2*N*b**3/(1+s)**3
def Phi3b(x,N): return 2*N*(1+3*x+3*x**2)/(x**3*(1+rt2*x+x*x)**arb("1.5"))
def Phi4b(x,N): return 6*N*(1+4*x+6*x**2+4*x**3)/(x**4*(1+rt2*x+x*x)**2)
def gF(x): return Hc(b*x)/x
def logg(x):   # log g(x) via lgamma; analytic since Re(bx/2)=x/(2rt2)>0 for x>0
    s = b*x
    return (sqpi/4).log() + (s/2).lgamma() - (s/2+arb("1.5")).lgamma() - x.log()
def gapF(x):
    s = b*x
    return b*((s/2).digamma() - (s/2+arb("1.5")).digamma())/2 - 1/x

def x0_root(N):   # point N (arb exact); bisection, ReP1 strictly decreasing in x
    lo = arb(1); hi = arb(10)**12
    assert arb(ReP1(lo,N).lower()) > 0 and arb(ReP1(hi,N).upper()) < 0
    for _ in range(1000):
        mid = (lo+hi)/2
        fm = ReP1(mid,N)
        if arb(fm.lower()) > 0: lo = mid
        elif arb(fm.upper()) < 0: hi = mid
        else: break      # undecidable at working precision; keep certified brackets
    assert arb(ReP1(lo,N).lower()) > 0 and arb(ReP1(hi,N).upper()) < 0
    return lo, hi   # f(lo)>0>f(hi) certified

def x0_enclosure(N):   # N arb (point or box): x0(N) in returned box for all N in box
    N1 = arb(N.lower()); N2 = arb(N.upper())
    lo1, hi1 = x0_root(N1)
    lo2, hi2 = x0_root(N2)
    # ReP1(x,N) = -a + N G(x); root increasing in N; lo1 has f>0 at N1 => >0 for all N>=N1 (G>0 there)
    x = iv(lo1, hi2)
    assert arb(ReP1(arb(lo1), N1).lower()) > 0 and arb(ReP1(arb(hi2), N2).upper()) < 0
    return x

def mom(j, cc):  # int_R |u|^j e^{-cc u^2} du = Gamma((j+1)/2)/cc^{(j+1)/2}
    return (arb(j+1)/2).gamma()/cc**((j+1)/2)

def certified(N_in, J=24, c_w=arb("1.0")):
    N = arb(N_in)
    th = (pi*N).sqrt()
    x0 = x0_enclosure(N)
    x0lo = arb(x0.lower()); x0hi = arb(x0.upper())
    ys = (pi/(2*N)).sqrt()
    y0 = 1/x0
    om = (N/rt2)*(ys**2 - y0**2/(1+rt2*y0+y0**2))
    assert om > 0
    Ph0 = PhiF(x0, N)
    Re0 = Ph0.real
    g0 = gF(x0)
    g0m = abs(g0)
    P0 = arb(g0m.upper())*arb(Re0.upper()).exp()
    c2pt = abs(ReP2(x0,N))/2
    c2 = arb(c2pt.lower())
    assert c2 > 0
    sig = 1/(2*c2).sqrt()
    w = (N.log()+c_w).sqrt()
    wsig = w*sig
    d = wsig if arb(wsig.lower()) <= arb((x0lo/3).lower()) else x0lo/3   # box-safe choice
    xl = x0lo - d; xr = x0hi + d
    assert xl > 0
    # curvature lower bounds on half windows via CELL COVERINGS (rigor fix vs step_33)
    cR = []; cL = []
    for j in range(50):
        cellR = iv(x0lo + (xr-x0lo)*arb(j)/50, x0lo + (xr-x0lo)*arb(j+1)/50)
        v = ReP2(cellR, N); assert arb(v.upper()) < 0
        cR.append((-arb(v.upper()))/2)
        cellL = iv(xl + (x0hi-xl)*arb(j)/50, xl + (x0hi-xl)*arb(j+1)/50)
        v = ReP2(cellL, N); assert arb(v.upper()) < 0
        cL.append((-arb(v.upper()))/2)
    c2R = arb(min(x.lower() for x in cR)); c2L = arb(min(x.lower() for x in cL))
    s0 = b*x0
    p2 = -N*(1+2*s0)/(x0*x0*(1+s0)**2)
    assert p2.real < 0
    assert -p2.real/2 >= c2
    Ph3c = Phi3F(x0,N); Ph3 = abs(Ph3c)
    GAPc = gapF(x0);  G1x = abs(GAPc)
    # W1: odd first-order x i sin(om u)
    W1 = P0*arb(om.upper())*(Ph3/6*mom(4,c2) + G1x*mom(2,c2))
    # W2/W3: subinterval loop with midpoint eval + certified reach
    W2 = arb(0); W3 = arb(0); maxReR = None
    logg0 = logg(x0)
    for side in [1,-1]:
        for j in range(J):
            u0 = d*arb(j)/J; u1 = d*arb(j+1)/J
            um = (u0+u1)/2
            h  = (u1-u0)/2
            uc = iv(u0,u1)
            xcell = x0 + side*uc
            xm    = x0 + side*um
            xL = arb(xcell.lower())
            assert xL > 0
            tL = xL/(2*rt2)
            M3c = arb(Phi3b(xcell,N).upper())
            M4c = arb(Phi4b(xcell,N).upper())
            G1c = arb(((arb("0.75")*(1/tL**2 + 1/tL) + 1/xL)).upper())
            G2c = arb(((arb("0.375")*(2/tL**3 + 1/tL**2) + 1/xL**2)).upper())
            su = side*um
            Rm = PhiF(xm,N) - Ph0 - I1*om*su - p2*su**2/2
            R4m = Rm - Ph3c*su**3/6
            gam = logg(xm) - logg0
            ga2m = gam - GAPc*su
            Rreach  = M3c*u1**2/2*h
            R4reach = M4c*u1**3/6*h
            Rb  = abs(Rm) + Rreach
            R4b = abs(R4m) + R4reach
            Gb  = abs(gam) + G1c*h
            Ga2b= abs(ga2m) + G2c*u1*h
            gw = (-c2*u0**2).exp()
            du = u1-u0
            V = Rb + Gb
            ReRup = arb((Rm.real + Rreach).upper())
            maxReR = ReRup if maxReR is None else max(maxReR, ReRup)
            gcell = Hr(xL/rt2)/xL          # >= |g| on cell (|Gamma(x+iy)|<=Gamma(x))
            qfac = (ReRup.exp()*(gcell/arb(g0m.lower()))).upper()
            q = arb(max(arb(1).mid(), arb(qfac).mid()))   # >= max(1, qfac)
            W2 += P0*gw*(R4b + Ga2b)*du
            W3 += P0*gw*(V*V/2)*q*du
    # model tail (IBP), both sides
    p2m = arb(abs(p2).lower()); omh = arb(om.upper())
    den = p2m*d - omh
    assert den > 0
    TM = P0*2*( (-c2*d*d).exp()/den + abs(p2)*(-c2*d*d).exp()/(den*den*2*c2*d) )
    # true tails
    GR = Hr(xr/rt2)/xr
    lamR = -ReP1(xr,N); assert lamR > 0
    TR = GR*arb(Re0.upper()).exp()*(-c2R*d*d).exp()/arb(lamR.lower())
    x2 = x0lo/2
    GLm2 = Hr(x2/rt2)/x2
    lamL = ReP1(xl,N); assert lamL > 0
    TLm = GLm2*arb(Re0.upper()).exp()*(-c2L*d*d).exp()/arb(lamL.lower())
    lam2 = ReP1(x2,N); assert lam2 > 0
    G1h = Hr(1/rt2)
    TLf = G1h*arb(RePhi(x2,N).upper()).exp()/arb(lam2.lower())
    cc0 = (2+rt2)/4
    TL0 = rt2*(1/(2*cc0.sqrt()))**N/(N-1)
    E_I = W1+W2+W3+TM+TR+TLm+TLf+TL0
    sc = N**arb("-1.5")*(-th).exp()
    MV = g0*Ph0.exp()*(-2*pi/p2).sqrt()*(om**2/(2*p2)).exp()
    Lt = -(I1*pi**arb("1.5")/2)*(-th+I1*th).exp()/N
    rho = abs(MV/Lt - 1)
    D = sqpi*rho*(-th).exp()/N
    C = ((2/pi)*E_I + D)/sc
    return dict(N=str(N_in), C=C, E_I=E_I, sc=sc, W1=W1, W2=W2, W3=W3, TM=TM,
                TR=TR, TLm=TLm, TLf=TLf, TL0=TL0, D=D, rho=rho, d=d, x0=x0,
                om=om, maxReR=maxReR, P0=P0, c2=c2, sig=sig, w=w)

if __name__ == "__main__":
    out = {}
    for Nv in ([int(x) for x in sys.argv[1:]] or [100, 400]):
        r = certified(Nv)
        out[str(Nv)] = {k: str(v) for k,v in r.items()}
        sc = r['sc']
        print(f"=== N={Nv}: certified C = {r['C']}")
        tot = 0
        for k in ['W1','W2','W3','TM','TR','TLm','TLf','TL0','D']:
            ratio = r[k]/sc
            print(f"   {k:4s}/sc = [{ratio}]")
        print(f"   ((2/pi)(W1+..+TL0)+D)/sc = C above; maxReR = {r['maxReR']}")
        print(f"   d={r['d']}, x0={r['x0']}, om={r['om']}")
    json.dump(out, open("cert_single_N.json","w"), indent=1)
