import sympy as sp, mpmath as mp, pickle
mp.mp.dps=70
a1,a2,a3,a4,a5=sp.symbols('a1 a2 a3 a4 a5')
G=sp.sympify([l.split(':',1)[1].strip() for l in open('R_factors.txt') if 'deg 22' in l][0])
C4=sp.sympify(open('C4.txt').read()); Q7=sp.sympify(open('Q7p.txt').read())
Q11=pickle.load(open('Q11_poly.pkl','rb'))
def Q11_mp(pt):
    tot=mp.mpf(0)
    for e,c in Q11.items():
        pw=mp.mpf(int(c.numerator))/int(c.denominator)
        for i in range(5):
            if e[i]: pw*=pt[i]**e[i]
        tot+=pw
    return tot
def a4_roots_exact(a1v,a2v):
    """exact isolation of positive real roots of G(a1,a2,a4), rational a1,a2. Returns isolating intervals."""
    Gp=sp.Poly(G.subs({a1:sp.Rational(a1v),a2:sp.Rational(a2v)}), a4, domain=sp.QQ)
    ivs=Gp.intervals()  # list of ((lo,hi), mult)
    return [(lo,hi) for ((lo,hi),m) in ivs if hi>0]
def refine_root(Gp, lo, hi, digits=60):
    """bisect to refine single root in (lo,hi) assuming sign change."""
    glo=sp.sign(Gp.subs(a4,lo)); 
    lo=mp.mpf(lo); hi=mp.mpf(hi)
    f=lambda x: mp.mpf(str(Gp.subs(a4, mp.nstr(x,50))))
    for _ in range(300):
        mid=(lo+hi)/2
        if mp.sign(f(mid))==mp.sign(f(lo)): lo=mid
        else: hi=mid
        if hi-lo<mp.mpf(10)**(-(digits+5)): break
    return (lo+hi)/2
def v10_hp(a1v,a2v):
    Gp=sp.Poly(G.subs({a1:sp.Rational(a1v),a2:sp.Rational(a2v)}), a4, domain=sp.QQ)
    out=[]
    C4s=C4.subs({a1:sp.Rational(a1v),a2:sp.Rational(a2v),a3:1})
    Q7s=Q7.subs({a1:sp.Rational(a1v),a2:sp.Rational(a2v),a3:1})
    for (lo,hi) in a4_roots_exact(a1v,a2v):
        if hi<=0: continue
        a4r=refine_root(Gp, lo, hi)
        if a4r<=0: continue
        # a5 common roots: C4 in a5 with a4=a4r
        C4p=sp.Poly(C4s.subs(a4, sp.Float(mp.nstr(a4r,60))), a5)
        cc4=[mp.mpf(str(c)) for c in C4p.all_coeffs()]
        try: r5=mp.polyroots(cc4, maxsteps=400)
        except Exception: 
            continue
        for s5 in r5:
            if abs(mp.im(s5))<mp.mpf('1e-40') and mp.re(s5)>mp.mpf('1e-55'):
                q7v=abs(Q7s.subs({a4:sp.Float(mp.nstr(a4r,50)),a5:sp.Float(mp.nstr(mp.re(s5),50))}))
                if q7v<mp.mpf('1e-35'):
                    out.append((a4r,mp.re(s5)))
    return out
for (x,y,lbl) in [((sp.Rational(1,100)), sp.Rational(1,100), "negative-region"),
                  ((sp.Rational(1,100)), sp.Rational(3), "positive-region")]:
    pts=v10_hp(x,y)
    print(f"(a1,a2)=({x},{y}) [{lbl}]: {len(pts)} V10+ points")
    for a4r,a5r in pts:
        q=Q11_mp([mp.mpf(x),mp.mpf(y),mp.mpf(1),a4r,a5r])
        print(f"  a4={mp.nstr(a4r,10)} a5={mp.nstr(a5r,10)} Q11={mp.nstr(q,10)} sign={'+' if q>0 else '-'}")
