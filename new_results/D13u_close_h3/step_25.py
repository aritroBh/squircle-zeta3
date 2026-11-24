25: def h3encl(T):  T = arb(T).upper();  return iv(rt2/3/(1+rt2*T+T*T), rt2/3)
59:     mS = (N*y0**3/2)*h3encl(y0)
60:     eRe0 = arb(mS.exp().upper())          # = exp(RePhi0+th) upper
70:     P0t = g0m_up*eRe0                     # -> W1, W2, W3, TM
134:    TR = GR*eRe0*.../arb(lamR.lower());  140: TLm = GLm2*eRe0*...
182:    z = (2*N*I1/pi**arb("1.5"))*g0box*mS.exp()*((-2*pi*invp2).sqrt())*((om**2*invp2/2).exp())*eiphi
183:    rhoz = arb(abs(z-1).upper());  184: D = sqpi*rhoz*...
197:    for k in range(nsub):  sig = iv(arb(k)/nsub, arb(k+1)/nsub);  cell(y, sig, ...)
