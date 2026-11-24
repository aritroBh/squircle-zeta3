import mpmath as mp
mp.mp.dps = 50
def P_direct_ab(a, b, p):
    e = 1/mp.mpf(p)
    f = lambda t: mp.sqrt((mp.mpf(a)*e*t**(e-1))**2 + (mp.mpf(b)*e*(1-t)**(e-1))**2)
    return 4*mp.quad(f, [mp.mpf(0), mp.mpf('0.5'), mp.mpf(1)])
def P_int(a, b, p, half):
    eta = 1/mp.mpf(p) - 1
    r = mp.sqrt(mp.mpf(a)**2+mp.mpf(b)**2); q = (mp.mpf(a)**2-mp.mpf(b)**2)/(mp.mpf(a)**2+mp.mpf(b)**2)
    ex = eta/2 if half else eta
    f = lambda t: (t*(1-t))**ex * mp.sqrt(mp.cosh(eta*mp.log(t/(1-t))) + q*mp.sinh(eta*mp.log(t/(1-t))))
    return 4*r*(1+eta)*mp.quad(f, [mp.mpf(0), mp.mpf('0.5'), mp.mpf(1)])
for (a,b,p) in [(1,1,'1.05'),(1,2,'1.05'),(1,3,'1.1'),(2,3,'0.97')]:
    pd = P_direct_ab(a,b,p)
    p_eta = P_int(a,b,p,False); p_half = P_int(a,b,p,True)
    print(f"(a,b,p)=({a},{b},{p}): direct={mp.nstr(pd,20)}")
    print(f"   exponent eta:   diff={mp.nstr(p_eta-pd,3)}")
    print(f"   exponent eta/2: diff={mp.nstr(p_half-pd,3)}")
