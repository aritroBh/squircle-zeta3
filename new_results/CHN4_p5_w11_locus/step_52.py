import sympy as sp, mpmath as mp, pickle
mp.mp.dps=60
a1,a2,a3,a4,a5=sp.symbols('a1 a2 a3 a4 a5')
G=sp.sympify([l.split(':',1)[1].strip() for l in open('R_factors.txt') if 'deg 22' in l][0])
C4=sp.sympify(open('C4.txt').read()); Q7=sp.sympify(open('Q7p.txt').read())
Q11=pickle.load(open('Q11_poly.pkl','rb'))
def to_mpf(r):  # exact rational -> mpf
    r=sp.Rational(r); return mp.mpf(int(r.p))/mp.mpf(int(r.q))
def Q11_mp(pt):
    tot=mp.mpf(0)
    for e,c in Q11.items():
        pw=mp.mpf(int(c.numerator))/int(c.denominator)
        for i in range(5):
            if e[i]: pw*=pt[i]**e[i]
        tot+=pw
    return tot
def v10_hp(a1v,a2v, a5min=mp.mpf('1e-40')):
    x,y=sp.Rational(a1v),sp.Rational(a2v)
    Gp=sp.Poly(G.subs({a1:x,a2:y}), a4, domain=sp.QQ)
    g_mp=sp.lambdify(a4, Gp.as_expr(), 'mpmath')
    C4s=C4.subs({a1:x,a2:y,a3:1}); Q7s=Q7.subs({a1:x,a2:y,a3:1})
    out=[]
    for (lo,hi),m in Gp.intervals():
        if to_mpf(hi)<=0: continue
        lo_m,hi_m=to_mpf(lo),to_mpf(hi)
        # refine by bisection (simple root) 
        flo=g_mp(lo_m)
        a,b=lo_m,hi_m
        for _ in range(400):
            mid=(a+b)/2
            if mp.sign(g_mp(mid))==mp.sign(flo): a=mid
            else: b=mid
            if b-a<mp.mpf(10)**(-55): break
        a4r=(a+b)/2
        if a4r<=0: continue
        # a5 common roots
        C4p=sp.Poly(C4s.subs(a4, sp.Float(mp.nstr(a4r,55))), a5)
        cc4=[mp.mpf(str(c)) for c in C4p.all_coeffs()]
        try: r5=mp.polyroots(cc4, maxsteps=400)
        except Exception: continue
        for s5 in r5:
            if abs(mp.im(s5))<mp.mpf('1e-40') and mp.re(s5)>a5min:
                q7v=Q7s.subs({a4:sp.Float(mp.nstr(a4r,50)),a5:sp.Float(mp.nstr(mp.re(s5),50))})
                if abs(q7v)<mp.mpf('1e-30'):
                    out.append((a4r,mp.re(s5)))
    return out
# test both regions
for (x,y,lbl) in [(sp.Rational(1,100),sp.Rational(1,100),"neg"),(sp.Rational(1,100),sp.Rational(3),"pos"),(sp.Rational(1),sp.Rational(1),"neg(1,1)")]:
    pts=v10_hp(x,y)
    print(f"(a1,a2)=({x},{y}) [{lbl}]: {len(pts)} pts")
    for a4r,a5r in pts:
        q=Q11_mp([mp.mpf(x),mp.mpf(y),mp.mpf(1),a4r,a5r])
        print(f"   a4={mp.nstr(a4r,10)} a5={mp.nstr(a5r,10)} Q11={mp.nstr(q,8)} {'+' if q>0 else '-'}")
