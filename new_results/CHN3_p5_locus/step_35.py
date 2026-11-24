exec(open('common.py').read())
R = sp.sympify(open('R_a1a2a4.txt').read())
G = sp.expand(R / (16*a1*(a1+1)))
C4,Q7 = get_C4_Q7()
g11 = G.subs({a1:1,a2:1})
nr = sp.nroots(g11, n=40, maxsteps=300)
posroots = sorted(set(sp.re(r) for r in nr if abs(sp.im(r))<sp.Float('1e-30') and sp.re(r)>0))
for r in posroots:
    if abs(r-1) < sp.Float('1e-20'): continue
    c4roots = sp.nroots(C4.subs({a1:1,a2:1,a3:1,a4:r}), n=40, maxsteps=300)
    q7roots = sp.nroots(Q7.subs({a1:1,a2:1,a3:1,a4:r}), n=40, maxsteps=300)
    # find matching pair
    best=None
    for x in c4roots:
        for y in q7roots:
            d=abs(x-y)
            if best is None or d<best[0]: best=(d,x,y)
    print(f"a4={sp.N(r,15)}: closest common a5 root diff={sp.N(best[0],3)}, a5={sp.N(best[1],25)}")
