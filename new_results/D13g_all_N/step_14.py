
src = open('step_37.py').read()
src = src.replace("# step_37.py - D13f:", "# step_38.py - D13g: step_37 + REGULAR TLf (kills N-box interval dependency) + tighter logcf form\n# D13f:")

old = """    eta2 = iv((x2/x0hi).lower(), (x2/x0lo).upper())
    lam2 = -a + (pi*rho/2)*Gh(eta2, rho, y)
    assert arb(lam2.lower()) > 0, "lam2 lost sign"
    G1h = arb(Hr(1/rt2).upper())
    x2i = iv(x2.lower(), x2.upper())
    vf = rt2/x2i + 1/x2i**2
    ReX2 = -a*x2i + th - (N/2)*(1+vf).log()      # = RePhi(x2)+th exactly
    TLf = G1h*arb(ReX2.upper()).exp()/arb(lam2.lower())"""
new = """    eta2 = iv((x2/x0hi).lower(), (x2/x0lo).upper())
    lam2 = -a + (pi*rho/2)*Gh(eta2, rho, y)
    assert arb(lam2.lower()) > 0, "lam2 lost sign"
    G1h = arb(Hr(1/rt2).upper())
    # REGULAR TLf: ReX2 = RePhi(x2)+th = th*(M1 - L2), both factors regular in (y,rho,eta2)
    #   M1 = 1 - eta2/(2 rho);  L2 = (rt2/4) ureg logcf(u), ureg = rt2 rho/eta2 + rho^2 y/eta2^2, u = y ureg
    M1 = 1 - eta2/(2*rho)
    ureg = rt2*rho/eta2 + rho*rho*y/(eta2*eta2)
    u = y*ureg
    logcf = ((1+u).log())/u                       # = log(1+u)/u, box-safe since u>0
    L2 = (rt2/4)*ureg*logcf
    MmL = M1 - L2
    assert arb(MmL.upper()) < 0, "TLf reg form lost sign (M1-L2 not < 0)"
    ReX2up = arb((arb(th.lower())*arb(MmL.upper())).upper())   # th>0, MmL<0
    TLf = G1h*arb(ReX2up.exp().upper())/arb(lam2.lower())"""
assert old in src
src = src.replace(old, new)
src = src.replace("def certified_y", "def certified_y38")
open('step_38.py','w').write(src)
print("step_38.py written", len(src))
