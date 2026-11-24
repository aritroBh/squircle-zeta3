import sympy as sp
x,y=sp.symbols('x y')
P=[1, None, 1-x,2-x,5-6*x+x**2,16-13*x+x**2,61-89*x+29*x**2-x**3,272-297*x+60*x**2-x**3,
   1385-2396*x+1134*x**2-124*x**3+x**4,7936-10841*x+3651*x**2-251*x**3+x**4,50521-101187*x+62372*x**2-12212*x**3+507*x**4-x**5]
P[1]=sp.Integer(1)  # sec+tan: A_1=1; check consistency at x=1: 1+sqrt2 tan(y/sqrt2) -> P1=1
N=11
for xv in [sp.Rational(1,2), sp.Rational(1,3)]:
    G=sum(sp.sympify(P[d]).subs(x,xv)*y**d/sp.factorial(d) for d in range(N))
    Gp=sp.diff(G,y)
    for (I,J) in [(0,2),(0,3),(1,2),(2,2),(0,4),(1,3)]:
        terms=[(i,j) for i in range(I+1) for j in range(J+1)]
        cs=sp.symbols('c0:%d'%len(terms))
        expr=sp.expand(Gp-sum(c*y**i*G**j for c,(i,j) in zip(cs,terms)))
        eqs=[expr.coeff(y,n) for n in range(N-1)]
        sol=sp.solve(eqs,cs,dict=True)
        print(xv,(I,J),len(terms),'eqs',N-1,'sol' if sol else 'none', sol[:1] if sol and len(terms)<N-2 else '')
