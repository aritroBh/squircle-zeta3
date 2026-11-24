exec(open('common.py').read())
R = sp.sympify(open('R_a1a2a4.txt').read())
S0 = sp.sympify(open('subres0.txt').read())
print("S0 == R:", sp.expand(S0-R)==0, " S0 == -R:", sp.expand(S0+R)==0)
A_ = sp.sympify(open('subres1_A.txt').read()); B_ = sp.sympify(open('subres1_B.txt').read())
C4,Q7 = get_C4_Q7()
# numeric check of a5=-B/A at second positive point a4 in (5,6)
g11 = G_expr = sp.expand(R/(16*a1*(a1+1))).subs({a1:1,a2:1})
nr = sp.nroots(g11, n=40, maxsteps=300)
for r in nr:
    if abs(sp.im(r))<sp.Float('1e-30') and sp.re(r)>5:
        a4v = sp.re(r)
        a5v = -B_.subs({a1:1,a2:1,a4:a4v})/A_.subs({a1:1,a2:1,a4:a4v})
        print(f"a4={sp.N(a4v,12)}: -B/A = {sp.N(a5v,20)}")
        # residual check: C4,Q7 at (1,1,1,a4v,a5v)
        print("  C4 resid:", sp.N(C4.subs({a1:1,a2:1,a3:1,a4:a4v,a5:a5v}),3),
              " Q7 resid:", sp.N(Q7.subs({a1:1,a2:1,a3:1,a4:a4v,a5:a5v}),3))
