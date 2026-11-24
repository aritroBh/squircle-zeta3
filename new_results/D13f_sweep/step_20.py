
import os
os.chdir("research/08_theory/D13f_sweep")
code = r'''# step_37.py - D13f: FIXED certified box assembly for
#   |M_k^main + (sqrt(pi)/N) e^{-th} cos th| <= C N^{-3/2} e^{-th}
# Fixes vs step_36 (the D13e-identified gap):
#  * gamma point/box evals (complex digamma/lgamma blow up over boxes in arb) replaced by the
#    PROVED real bounds |g'|<= (3/4)(1/t^2+1/t)+1/x, |g''|<= (3/8)(2/t^3+1/t^2)+1/x^2, t=x/(2 sqrt2)
#  * lambdaR/lambdaL/lambda2 by the PROVED regular form  -ReP1(eta x0,N) = a - (pi rho/2) Gh(eta)
#  * g0 (complex) via midpoint reference + propagation |g(x)/g(xref)| <= exp(G1 |x-xref|)
#  * |Phi3|<=Phi3b, |Phi4|<=Phi4b (PROVED monotone bounds) instead of cancelling complex evals
#  * q-factor uses |g| DECREASING on x>0 (PROVED: Re g'(x) < 0 termwise)
from flint import arb, acb, ctx
import json, sys, time
ctx.dps = 70
pi = arb.pi(); rt2 = arb(2).sqrt(); a = pi/(2*rt2); sqpi = pi.sqrt()
b = (acb(0,1)*pi/4).exp(); I1 = acb(0,1)
def iv(lo, hi):
    lo = arb(lo).lower(); hi = arb(hi).upper()
    return arb((lo+hi)/2, (hi-lo)/2)
def Hr(t): return sqpi/4*(t/2).gamma()/(t/2+arb("1.5")).gamma()
def Hc(s): return sqpi/4*(s/2).gamma()/(s/2+arb("1.5")).gamma()
def gF(x): return Hc(b*x)/x
def Phi3b(x,N): return 2*N*(1+3*x+3*x**2)/(x**3*(1+rt2*x+x*x)**arb("1.5"))
def Phi4b(x,N): return 6*N*(1+4*x+6*x**2+4*x**3)/(x**4*(1+rt2*x+x*x)**2)
def mom(j, cc): return (arb(j+1)/2).gamma()/cc**((arb(j)+1)/2)
def h3encl(T):
    T = arb(T).upper()
    return iv(rt2/3/(1+rt2*T+T*T), rt2/3)
def curv2(y, rho, eta):   # = |RePhi''(eta x0)|/2 exactly
    num = rt2*eta**3*rho + 4*eta**2*rho**2*y + 2*rt2*eta*rho**3*y**2 + rho**4*y**3
    den = eta**2*(eta**2+rt2*eta*rho*y+rho**2*y**2)**2
    return (pi*rho**2*y/4)*num/den
def Gh(eta, rho, y):      # -ReP1(eta*x0,N) = a - (pi rho/2) Gh(eta)
    t = rho*y
    return rho*(eta+t*rt2)/(rt2*eta*(eta*eta+rt2*eta*t+t*t))
def G1bound(x):           # >= |gamma'(x)|, decreasing in x
    t = x/(2*rt2)
    return arb("0.75")*(1/t**2 + 1/t) + 1/x
def G2bound(x):           # >= |gamma''(x)|, decreasing in x
    t = x/(2*rt2)
    return arb("0.375")*(2/t**3 + 1/t**2) + 1/x**2

def cell(y, sig, J, c_w, coll):
    rho = 1 + y*y*sig
    y0 = rho*y
    x0 = 1/y0
    x0lo = arb(x0.lower()); x0hi = arb(x0.upper())
    N = pi/(2*y*y); th = pi/(rt2*y)
    den1 = 1+rt2*rho*y+rho*rho*y*y
    Omn = rt2*rho + rho*rho*y - 2*y*sig - y**3*sig**2
    om = (pi/(2*rt2))*y*Omn/den1
    assert om > 0, "om lost sign"
    mS = (N*y0**3/2)*h3encl(y0)
    eRe0 = arb(mS.exp().upper())          # = exp(RePhi0+th) upper
    xref = arb(x0.mid())
    hx = (x0hi-x0lo)/2
    G1box = arb(G1bound(x0lo).upper())
    rg = G1box*hx
    g0ref = gF(xref)
    g0box = g0ref * acb(arb(0,rg), arb(0,rg)).exp()
    g0m_up_prop = arb(abs(g0ref).upper())*arb(rg.exp().upper())
    g0m_up_beta = arb((Hr(x0lo/rt2)/x0lo).upper())
    g0m_up = min(g0m_up_prop, g0m_up_beta)
    P0t = g0m_up*eRe0
    t = rho*y
    numR = t**3 + 2*rt2*t**2 + 4*t + rt2
    D2 = t**2 + 2*rt2*t + 4
    p2 = -(pi/2)*rho**3*y*acb(t+rt2, rt2)/acb(t**2+rt2*t, rt2*t+1)
    Rep2 = -(pi/2)*rho**3*y*numR/(1+rt2*t+t**2)**2
    assert Rep2 < 0
    c2  = arb((-Rep2/2).lower()); assert c2 > 0
    c2up= arb((-Rep2/2).upper())
    p2mag = (pi/2)*rho**3*y*(D2.sqrt())/(t**2+rt2*t+1)
    invp2 = -(2/(pi*rho**3*y))*acb(numR, t+rt2)/D2
    sg = 1/(2*c2).sqrt()
    w = (N.log()+c_w).sqrt()
    wsig = w*sg
    d = wsig if arb(wsig.upper()) <= arb((x0lo/3).lower()) else x0lo/3
    dlo = arb(d.lower()); dup = arb(d.upper())
    xl = x0lo - d; xr = x0hi + d
    assert arb(xl.lower()) > 0
    etalo = arb((xl/x0hi).lower()); etahi = arb((xr/x0lo).upper())
    cR = []; cL = []
    for j in range(60):
        eR = iv(arb(1)+(etahi-1)*arb(j)/60, arb(1)+(etahi-1)*arb(j+1)/60)
        cR.append(arb(curv2(y, rho, eR).lower()))
        eL = iv(etalo+(1-etalo)*arb(j)/60, etalo+(1-etalo)*arb(j+1)/60)
        cL.append(arb(curv2(y, rho, eL).lower()))
    c2R = min(cR); c2L = min(cL)
    assert c2R > 0 and c2L > 0
    Ph3B = arb(Phi3b(x0lo, arb(N.upper())).upper())
    G1B  = arb(G1bound(x0lo).upper())
    W1 = P0t*arb(om.upper())*(Ph3B/6*mom(4,c2) + G1B*mom(2,c2))
    W2 = arb(0); W3 = arb(0)
    for side in [1,-1]:
        c2side = c2R if side==1 else c2L
        for j in range(J):
            u0 = d*arb(j)/J; u1 = d*arb(j+1)/J
            u1up = arb(u1.upper()); u0lo = arb(u0.lower())
            xsegL = x0lo if side==1 else x0lo - u1up
            assert xsegL > 0
            M3 = arb(Phi3b(xsegL, arb(N.upper())).upper())
            M4 = arb(Phi4b(xsegL, arb(N.upper())).upper())
            G1 = arb(G1bound(xsegL).upper())
            G2 = arb(G2bound(xsegL).upper())
            R4b = M4*u1up**4/24
            Ga2b= G2*u1up**2/2
            Rb  = M3*u1up**3/6
            Gb  = G1*u1up
            V = Rb + Gb
            qe = (c2up-c2side)*u1up**2
            if side==-1: qe = qe + G1*u1up
            q = arb(max(arb(1).mid(), arb(qe.exp().upper())))
            gw = arb((-c2*u0lo**2).exp().upper())
            du = arb((u1-u0).upper())
            W2 += P0t*gw*(R4b+Ga2b)*du
            W3 += P0t*gw*(V*V/2)*q*du
    p2m = arb(p2mag.lower()); omh = arb(om.upper())
    den = p2m*dlo - omh
    assert arb(den.lower()) > 0, "TM den lost sign"
    ecd = arb((-c2*dlo**2).exp().upper())
    TM = P0t*2*( ecd/den + arb(p2mag.upper())*ecd/(arb(den.lower())**2*2*c2*dlo) )
    xrt = x0lo + dlo
    GR = arb((Hr(xrt/rt2)/xrt).upper())
    etaR = iv((xr/x0hi).lower(), (xr/x0lo).upper())
    lamR = a - (pi*rho/2)*Gh(etaR, rho, y)
    assert arb(lamR.lower()) > 0, "lamR lost sign"
    TR = GR*eRe0*arb((-c2R*dlo**2).exp().upper())/arb(lamR.lower())
    x2 = x0lo/2
    GLm2 = arb((Hr(x2/rt2)/x2).upper())
    etaL = iv((xl/x0hi).lower(), (xl/x0lo).upper())
    lamL = -a + (pi*rho/2)*Gh(etaL, rho, y)
    assert arb(lamL.lower()) > 0, "lamL lost sign"
    TLm = GLm2*eRe0*arb((-c2L*dlo**2).exp().upper())/arb(lamL.lower())
    eta2 = iv((x2/x0hi).lower(), (x2/x0lo).upper())
    lam2 = -a + (pi*rho/2)*Gh(eta2, rho, y)
    assert arb(lam2.lower()) > 0, "lam2 lost sign"
    G1h = arb(Hr(1/rt2).upper())
    x2i = iv(x2.lower(), x2.upper())
    vf = rt2/x2i + 1/x2i**2
    ReX2 = -a*x2i + th - (N/2)*(1+vf).log()      # = RePhi(x2)+th exactly
    TLf = G1h*arb(ReX2.upper()).exp()/arb(lam2.lower())
    cc0 = (2+rt2)/4
    ex0 = (1/(2*cc0.sqrt())).log() + (pi/N).sqrt()
    assert arb(ex0.upper()) < 0
    TL0 = rt2*arb((ex0*N).upper()).exp()/(N-1)
    E_I = W1+W2+W3+TM+TR+TLm+TLf+TL0
    zz = rt2*rho*y/(2+rt2*rho*y)
    z2 = zz*zz
    A1m = iv(0, (z2/(1+z2)).upper())
    T2n = 2*y*y*sig - rt2*rho*y - 2*rho*A1m*y
    Bh = -y*sig/(2*rho) + T2n/(2*y*(2+rt2*rho*y))
    phireg = (pi/rt2)*Bh
    eiphi = (I1*phireg).exp()
    z = (2*N*I1/pi**arb("1.5"))*g0box*mS.exp()*((-2*pi*invp2).sqrt())*((om**2*invp2/2).exp())*eiphi
    rhoz = arb(abs(z-1).upper())
    D = sqpi*rhoz*arb(N.sqrt().upper())
    N15 = arb(N**arb("1.5").upper())
    C = (2/pi)*E_I*N15 + D
    coll['C'] = max(coll['C'], arb(C.upper()))
    for k,v in [('W1',W1),('W2',W2),('W3',W3),('TM',TM),('TR',TR),('TLm',TLm),('TLf',TLf),('TL0',TL0),('D',D)]:
        coll[k] = max(coll[k], arb((v*N15).upper()) if k!='D' else arb(v.upper()))
    return arb(C.upper())

def certified_y(y1, y2, J=24, nsub=1, c_w=arb("1.0")):
    y = iv(y1, y2)
    assert arb(y2) <= arb("0.78"), "y0 enclosure proved only for y<=0.78"
    coll = dict(C=arb(0), W1=arb(0), W2=arb(0), W3=arb(0), TM=arb(0), TR=arb(0),
                TLm=arb(0), TLf=arb(0), TL0=arb(0), D=arb(0))
    for k in range(nsub):
        sig = iv(arb(k)/nsub, arb(k+1)/nsub)
        cell(y, sig, J, c_w, coll)
    return coll

if __name__ == "__main__":
    import math
    ys = lambda N: math.sqrt(math.pi/(2*N))
    out = {}
    for Nv in [100, 400]:
        yv = ys(Nv)
        t0=time.time()
        r = certified_y(yv*(1-1e-6), yv*(1+1e-6), nsub=4)
        print(f"N~{Nv}: C = {r['C']}   ({time.time()-t0:.1f}s)")
        print("   pieces*N^1.5:", {k: round(float(v.mid()),4) for k,v in r.items() if k!='C'})
    cells = [(100,150),(150,225),(225,340),(340,400)]
    for (n1,n2) in cells:
        t0=time.time()
        r = certified_y(ys(n2), ys(n1), nsub=4)
        print(f"N-box [{n1},{n2}]: C = {r['C'].upper()}  ({time.time()-t0:.1f}s)")
        print("   pieces*N^1.5:", {k: round(float(v.mid()),3) for k,v in r.items() if k!='C'})
        out[f"{n1}-{n2}"] = {k: str(v) for k,v in r.items()}
        json.dump(out, open("cert_boxes37.json","w"), indent=1)
'''
with open("step_37.py","w") as f:
    f.write(code)
import subprocess
PY = "<home>)/.venv/bin/python"
r = subprocess.run([PY, "-c", "import ast; ast.parse(open('step_37.py').read()); print('syntax ok')"],
                   capture_output=True, text=True)
print(r.stdout, r.stderr[:1000])
