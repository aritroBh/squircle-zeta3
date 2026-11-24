import sympy as sp, mpmath as mp, pickle
mp.mp.dps=80
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
def v10_points_hp(a1v,a2v):
    """high-precision V10+ points at rational (a1,a2). Returns list of (a4,a5) mpf."""
    Gp=sp.Poly(G.subs({a1:sp.Rational(a1v),a2:sp.Rational(a2v)}), a4)
    cc=[mp.mpf(int(c)) for c in Gp.all_coeffs()]
    try: roots=mp.polyroots(cc, maxsteps=200)
    except Exception as ex:
        # fallback: sympy nroots
        nr=sp.nroots(Gp, n=80, maxsteps=300)
        roots=[mp.mpf(str(sp.re(r)))+1j*mp.mpf(str(sp.im(r))) for r in nr]
    out=[]
    C4s=C4.subs({a1:sp.Rational(a1v),a2:sp.Rational(a2v),a3:1})
    Q7s=Q7.subs({a1:sp.Rational(a1v),a2:sp.Rational(a2v),a3:1})
    for r in roots:
        if abs(mp.im(r))<mp.mpf('1e-40') and mp.re(r)>0:
            a4r=mp.re(r)
            # common a5: roots of C4 in a5
            C4p=sp.Poly(C4s.subs(a4, mp.nstr(a4r,60)), a5)
            cc4=[mp.mpf(str(c)) for c in C4p.all_coeffs()]
            try: r5=mp.polyroots(cc4, maxsteps=200)
            except Exception: continue
            for s5 in r5:
                if abs(mp.im(s5))<mp.mpf('1e-30') and mp.re(s5)>mp.mpf('1e-50'):
                    q7val=abs(Q7s.subs({a4:mp.nstr(a4r,50),a5:mp.nstr(mp.re(s5),50)}))
                    if q7val<mp.mpf('1e-30'):
                        out.append((a4r,mp.re(s5)))
    return out
for (x,y,lbl) in [((sp.Rational(1,100)), sp.Rational(3), "positive-region"),
                  ((sp.Rational(1,100)), sp.Rational(1,100), "negative-region")]:
    pts=v10_points_hp(x,y)
    print(f"\n(a1,a2)=({x},{y}) [{lbl}]: {len(pts)} V10+ points")
    for a4r,a5r in pts:
        q=Q11_mp([mp.mpf(x),mp.mpf(y),mp.mpf(1),a4r,a5r])
        print(f"  a4={mp.nstr(a4r,12)} a5={mp.nstr(a5r,12)} Q11={mp.nstr(q,12)}  sign={'+' if q>0 else '-'}")
