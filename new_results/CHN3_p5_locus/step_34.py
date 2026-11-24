exec(open('common.py').read())
R = sp.sympify(open('R_a1a2a4.txt').read())
G = sp.expand(R / (16*a1*(a1+1)))
g11 = G.subs({a1:1,a2:1})
nr = sp.nroots(g11, n=30, maxsteps=200)
posroots = [r for r in nr if abs(sp.im(r))<sp.Float('1e-20') and sp.re(r)>0]
print("positive real roots a4 of G(1,1,a4):")
for r in posroots: print("  ", sp.N(r, 25))
C4,Q7 = get_C4_Q7()
# For each, find common a5 root of C4(1,1,1,a4,a5), Q7'
for r in posroots:
    a4v = sp.re(r)
    c4v = sp.Poly(C4.subs({a1:1,a2:1,a3:1,a4:a4v}), a5)
    q7v = sp.Poly(Q7.subs({a1:1,a2:1,a3:1,a4:a4v}), a5)
    g = sp.gcd(c4v, q7v, a5)
    print(f"a4={sp.N(a4v,12)}: gcd degree {g.degree()} , roots {sp.nroots(g, n=20) if g.degree()>0 else '-'}")
