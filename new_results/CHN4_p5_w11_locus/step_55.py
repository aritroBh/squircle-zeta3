import sympy as sp, mpmath as mp, pickle
mp.mp.dps=60
a1,a2,a3,a4,a5=sp.symbols('a1 a2 a3 a4 a5')
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
# Newton continuation in (a4,a5) for fixed (a1,a2)
def make_funcs(a1v,a2v):
    C4s=C4.subs({a1:sp.Rational(a1v),a2:sp.Rational(a2v),a3:1})
    Q7s=Q7.subs({a1:sp.Rational(a1v),a2:sp.Rational(a2v),a3:1})
    f1=sp.lambdify((a4,a5), C4s, 'mpmath'); f2=sp.lambdify((a4,a5), Q7s, 'mpmath')
    return f1,f2
def trace(a1v, a2_start, a2_end, a4_guess, a5_guess, nsteps=400):
    f1,f2=make_funcs(a1v,a2_start)
    # initial Newton at start
    sol=mp.findroot([f1,f2], (a4_guess,a5_guess), tol=mp.mpf('1e-45'), maxsteps=100)
    pts=[]
    da=(mp.mpf(a2_end)-mp.mpf(a2_start))/nsteps
    cur=(sol[0],sol[1])
    for i in range(nsteps+1):
        a2v=mp.mpf(a2_start)+i*da
        f1,f2=make_funcs(a1v,sp.Float(mp.nstr(a2v,50)))
        try:
            sol=mp.findroot([f1,f2], cur, tol=mp.mpf('1e-45'), maxsteps=100)
            cur=(sol[0],sol[1])
        except Exception as ex:
            print(f"  continuation failed at a2={mp.nstr(a2v,8)}: {ex}"); break
        q=Q11_mp([mp.mpf(sp.Rational(a1v)),a2v,mp.mpf(1),cur[0],cur[1]])
        pts.append((a2v,cur[0],cur[1],q))
    return pts
# start from the negative point (a1=1/100, a2=1/100, a4=0.2239, a5=0.9753), increase a2 to 3
pts=trace(sp.Rational(1,100), sp.Rational(1,100), sp.Rational(3), 
          mp.mpf('0.2239'), mp.mpf('0.9753'), nsteps=300)
print(f"traced {len(pts)} points")
# find sign change
prev=None
for (a2v,a4v,a5v,q) in pts:
    if prev is not None and prev[3]*q<0:
        print(f"SIGN CHANGE between a2={mp.nstr(prev[0],10)} (Q11={mp.nstr(prev[3],6)}) and a2={mp.nstr(a2v,10)} (Q11={mp.nstr(q,6)})")
        print(f"   bracket a4~{mp.nstr(a4v,8)} a5~{mp.nstr(a5v,8)}")
    prev=(a2v,a4v,a5v,q)
print("first/last:", [(mp.nstr(p[0],6), mp.nstr(p[3],6)) for p in [pts[0],pts[-1]]])
