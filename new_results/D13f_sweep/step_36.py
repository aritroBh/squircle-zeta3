from flint import arb, acb, ctx
import json, sys, time
ctx.dps = 70
pi = arb.pi(); rt2 = arb(2).sqrt(); a = pi/(2*rt2); sqpi = pi.sqrt()
b = (acb(0,1)*pi/4).exp()
I1 = acb(0,1)
def iv(lo, hi):
    lo = arb(lo).lower(); hi = arb(hi).upper()
    return arb((lo+hi)/2, (hi-lo)/2)
def Hr(t): return sqpi/4*(t/2).gamma()/(t/2+arb("1.5")).gamma()
def Hc(s): return sqpi/4*(s/2).gamma()/(s/2+arb("1.5")).gamma()
def PhiF(x, N):
    s = b*x
    return I1*pi*s/2 + N*(s/(1+s)).log()
def ReP1(x,N):  return -a + N/x - N*(x+1/rt2)/(1+rt2*x+x*x)
def Phi3F(x,N):
    s = b*x
    return 2*N/x**3 - 2*N*b**3/(1+s)**3
def Phi3b(x,N): return 2*N*(1+3*x+3*x**2)/(x**3*(1+rt2*x+x*x)**arb("1.5"))
def Phi4b(x,N): return 6*N*(1+4*x+6*x**2+4*x**3)/(x**4*(1+rt2*x+x*x)**2)
def gF(x): return Hc(b*x)/x
def logg(x):
    s = b*x
    return (sqpi/4).log() + (s/2).lgamma() - (s/2+arb("1.5")).lgamma() - x.log()
def gapF(x):
    s = b*x
    return b*((s/2).digamma() - (s/2+arb("1.5")).digamma())/2 - 1/x
def mom(j, cc): return (arb(j+1)/2).gamma()/cc**((j+1)/2)
def h3encl(T):
    T = arb(T).upper()
    return iv(rt2/3/(1+rt2*T+T*T), rt2/3)
# ReP2(x) = -N y0^2 A2e(eta), with -A2e = y*(rt2 eta^3 rho + 4 eta^2 rho^2 y + 2 rt2 eta rho^3 y^2 + rho^4 y^3)/(eta^2 (eta^2+rt2 eta rho y+rho^2 y^2)^2)
# => |ReP2|/2 = (pi rho^2 /4) * (rt2 eta^3 rho + 4 eta^2 rho^2 y + 2 rt2 eta rho^3 y^2 + rho^4 y^3)/(eta^2 (eta^2+rt2 eta rho y+rho^2 y^2)^2)
def curv2(y, rho, eta):
    num = rt2*eta**3*rho + 4*eta**2*rho**2*y + 2*rt2*eta*rho**3*y**2 + rho**4*y**3
    den = eta**2*(eta**2+rt2*eta*rho*y+rho**2*y**2)**2
    return (pi*rho**2*y/4)*num/den

def certified_y(y1, y2, J=24, c_w=arb("1.0")):
    y = iv(y1, y2)
    assert arb(y2) <= arb("0.78"), "enclosure y0 in [y,y(1+y^2)] proved only for y<=0.78 (N>=2.6)"
    sig = iv(0, 1)      # y0 = rho y, rho = 1 + y^2 sig: PROVED (F increasing, F(y)<=y^2, F(y(1+y^2))>=y^2 by polynomial ineq.)
    rho = 1 + y*y*sig
    y0 = rho*y
    x0 = 1/y0
    N = pi/(2*y*y)
    th = pi/(rt2*y)
    den1 = 1+rt2*rho*y+rho*rho*y*y
    Omn = rt2*rho + rho*rho*y - 2*y*sig - y**3*sig**2
    Om = y*Omn/den1
    om = (pi/(2*rt2))*Om
    assert om > 0
    mS = (N*y0**3/2)*h3encl(y0)
    eRe0 = mS.exp()
    g0 = gF(x0)
    g0m = abs(g0)
    P0t = arb(g0m.upper())*arb(eRe0.upper())
    # c2 at eta=1 (regular, no cancellation)
    c2i = curv2(y, rho, iv(1,1))
    t = rho*y
    numR = t**3 + 2*rt2*t**2 + 4*t + rt2          # Re/Im parts of p2, 1/p2 (cancellation-free)
    D2 = t**2 + 2*rt2*t + 4
    p2 = -(pi/2)*rho**3*y*acb(t+rt2, rt2)/acb(t**2+rt2*t, rt2*t+1)
    Rep2 = -(pi/2)*rho**3*y*numR/(1+rt2*t+t**2)**2
    assert Rep2 < 0
    c2 = arb((-Rep2/2).lower()); assert c2 > 0
    p2mag = (pi/2)*rho**3*y*(D2.sqrt())/(t**2+rt2*t+1)
    invp2 = -(2/(pi*rho**3*y))*acb(numR, t+rt2)/D2
    sg = 1/(2*c2).sqrt()
    w = (N.log()+c_w).sqrt()
    wsig = w*sg
    x0lo = arb(x0.lower()); x0hi = arb(x0.upper())
    d = wsig if arb(wsig.lower()) <= arb((x0lo/3).lower()) else x0lo/3
    xl = x0lo - d; xr = x0hi + d
    assert xl > 0
    # curvature lower bounds on half windows via eta cells (regular curv2)
    etalo = arb((xl/x0hi).lower()); etahi = arb((xr/x0lo).upper())
    cR = []; cL = []
    for j in range(50):
        eR = iv(arb(1)+(etahi-1)*arb(j)/50, arb(1)+(etahi-1)*arb(j+1)/50)
        cR.append(curv2(y, rho, eR))
        eL = iv(etalo+(1-etalo)*arb(j)/50, etalo+(1-etalo)*arb(j+1)/50)
        cL.append(curv2(y, rho, eL))
    c2R = arb(min(t.lower() for t in cR)); c2L = arb(min(t.lower() for t in cL))
    Ph3c = Phi3F(x0,N); Ph3 = abs(Ph3c)
    GAPc = gapF(x0);  G1x = abs(GAPc)
    W1 = P0t*arb(om.upper())*(Ph3/6*mom(4,c2) + G1x*mom(2,c2))
    W2 = arb(0); W3 = arb(0); maxReR = None
    logg0 = logg(x0)
    for side in [1,-1]:
        for j in range(J):
            u0 = d*arb(j)/J; u1 = d*arb(j+1)/J
            um = (u0+u1)/2; h = (u1-u0)/2
            uc = iv(u0,u1)
            xcell = x0 + side*uc
            xm    = x0 + side*um
            xL = arb(xcell.lower()); assert xL > 0
            tL = xL/(2*rt2)
            M3c = arb(Phi3b(xcell,N).upper())
            M4c = arb(Phi4b(xcell,N).upper())
            G1c = arb((arb("0.75")*(1/tL**2 + 1/tL) + 1/xL).upper())
            G2c = arb((arb("0.375")*(2/tL**3 + 1/tL**2) + 1/xL**2).upper())
            su = side*um
            Rm = PhiF(xm,N) - PhiF(x0,N) - I1*om*su - p2*su**2/2
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
            gcell = Hr(xL/rt2)/xL
            qfac = (ReRup.exp()*(gcell/arb(g0m.lower()))).upper()
            q = arb(max(arb(1).mid(), arb(qfac).mid()))
            W2 += P0t*gw*(R4b + Ga2b)*du
            W3 += P0t*gw*(V*V/2)*q*du
    p2m = arb(p2mag.lower()); omh = arb(om.upper())
    den = p2m*d - omh
    assert den > 0
    TM = P0t*2*( (-c2*d*d).exp()/den + p2mag*(-c2*d*d).exp()/(den*den*2*c2*d) )
    GR = Hr(xr/rt2)/xr
    lamR = -ReP1(xr,N); assert lamR > 0
    TR = GR*arb(eRe0.upper())*(-c2R*d*d).exp()/arb(lamR.lower())
    x2 = x0lo/2
    GLm2 = Hr(x2/rt2)/x2
    lamL = ReP1(xl,N); assert lamL > 0
    TLm = GLm2*arb(eRe0.upper())*(-c2L*d*d).exp()/arb(lamL.lower())
    lam2 = ReP1(x2,N); assert lam2 > 0
    G1h = Hr(1/rt2)
    x2i = iv(x2.lower(), x2.upper())
    vf = rt2/x2i + 1/x2i**2
    ReX2 = -a*x2i + th - (N/2)*(1+vf).log()
    TLf = G1h*arb(ReX2.upper()).exp()/arb(lam2.lower())
    cc0 = (2+rt2)/4
    ex0 = (1/(2*cc0.sqrt())).log() + (pi/N).sqrt()
    assert ex0 < 0
    TL0 = rt2*(ex0*N).exp()*N**arb("1.5")/(N-1)
    E_I = W1+W2+W3+TM+TR+TLm+TLf+TL0
    # D-term, sigma-correlated phase: phireg = (pi/rt2) * Bh
    zz = rt2*rho*y/(2+rt2*rho*y)
    z2 = zz*zz
    A1m = iv(0, (z2/(1+z2)).upper())       # 1 - atan(z)/z in [0, z^2/(1+z^2)]
    T2n = 2*y*y*sig - rt2*rho*y - 2*rho*A1m*y   # = 2 rho A1 - 2 - rt2 rho y, as y*(...)
    Bh = -y*sig/(2*rho) + T2n/(2*y*(2+rt2*rho*y))
    phireg = (pi/rt2)*Bh
    eiphi = (I1*phireg).exp()
    z = (2*N*I1/pi**arb("1.5"))*g0*mS.exp()*((-2*pi*invp2).sqrt())*((om**2*invp2/2).exp())*eiphi
    rhoz = abs(z-1)
    D = sqpi*rhoz*N.sqrt()
    C = (2/pi)*E_I*N**arb("1.5") + D
    return dict(ybox=(y1,y2), Nlo=str(N.lower()), Nhi=str(N.upper()), C=C,
        W1=W1*N**arb("1.5"), W2=W2*N**arb("1.5"), W3=W3*N**arb("1.5"), TM=TM*N**arb("1.5"),
        TR=TR*N**arb("1.5"), TLm=TLm*N**arb("1.5"), TLf=TLf*N**arb("1.5"), TL0=TL0, D=D,
        x0=x0, d=d, om=om, maxReR=maxReR)

if __name__ == "__main__":
    import math
    out = {}
    ys = lambda N: math.sqrt(math.pi/(2*N))
    for Nv in [100, 400]:
        yv = ys(Nv)
        r = certified_y(yv*(1-1e-6), yv*(1+1e-6))
        print(f"N~{Nv}: C = {r['C']}")
        print("   pieces/sc:", {k: round(float(v.mid()),4) for k,v in r.items() if k in ['W1','W2','W3','TM','TR','TLm','TLf','TL0','D']})
    cells = [(100,150),(150,225),(225,340),(340,400)]
    for (n1,n2) in cells:
        t0=time.time()
        r = certified_y(ys(n2), ys(n1))
        print(f"N-box [{n1},{n2}]: Cupper = {r['C'].upper()}  ({time.time()-t0:.1f}s)")
        print("   pieces/sc:", {k: round(float(v.mid()),3) for k,v in r.items() if k in ['W1','W2','W3','TM','TR','TLm','TLf','TL0','D']})
        out[f"{n1}-{n2}"] = {k: str(v) for k,v in r.items() if k not in ('ybox','Nlo','Nhi')}
        json.dump(out, open("cert_boxes.json","w"), indent=1)
