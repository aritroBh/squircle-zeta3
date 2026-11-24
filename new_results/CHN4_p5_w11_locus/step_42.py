import sympy as sp, time
a1,a2,a4=sp.symbols('a1 a2 a4')
G=sp.sympify([l.split(':',1)[1].strip() for l in open('R_factors.txt') if 'deg 22' in l][0])
G11=G.subs({a1:1,a2:1})
t0=time.time()
nr=sp.nroots(G11, n=30, maxsteps=200)
print(f"nroots time {time.time()-t0:.1f}s")
print("roots of G(1,1,a4):")
for r in nr:
    if abs(sp.im(r))<sp.Float('1e-20'):
        rv=sp.re(r)
        print("   ", rv, "POSITIVE" if rv>0 else "")
