exec(open('screen_R5.py').read().split("# (2) independent numeric")[0].split("x=sp.Symbol")[0])
src=open('screen_R5.py').read()
pre=src.split("tested=0;bad=[]")[0].split("# (2) independent numeric: gaussian-rational series mod s^6")[1]
exec(pre)
for dd in range(3,8):
    for pp in range((dd+1)//2,dd+1):
        ser=Sser(dd,pp)
        if ser[0].a==0 and ser[0].b==0: print(dd,pp,'S(t*)=0'); continue
        r=Rl(ser,5)
        pi=Pi(5,dd)
        if pi==0: print(dd,pp,'Pi5=0'); continue
        v=sp.expand(N5.subs({d_:dd,m_:pp})/(3*pi))
        print(dd,pp,'match' if (sp.Rational(sp.re(v))==r.a and sp.Rational(sp.im(v))==r.b) else 'MISMATCH')
