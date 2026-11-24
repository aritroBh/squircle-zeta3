import sympy as sp
x=sp.symbols('x')
P=[1, 1, 1-x,2-x,5-6*x+x**2,16-13*x+x**2,61-89*x+29*x**2-x**3,272-297*x+60*x**2-x**3,
   1385-2396*x+1134*x**2-124*x**3+x**4,7936-10841*x+3651*x**2-251*x**3+x**4,50521-101187*x+62372*x**2-12212*x**3+507*x**4-x**5]
# S-fraction (Stieltjes) via quotient-difference: f = c0/(1 - a1 y/(1 - a2 y/(...)))
def qd(c,n):
    # returns a_k via the qd algorithm on moments c
    e=[[sp.Integer(0)]*(n+2) for _ in range(n+2)]
    q=[[None]*(n+2) for _ in range(n+2)]
    for j in range(n): q[1][j]=sp.cancel(c[j+1]/c[j])
    out=[]
    for k in range(1,n):
        for j in range(n-2*k+1):
            e[k][j]=sp.cancel(q[k][j+1]-q[k][j]+e[k-1][j+1])
        for j in range(n-2*k):
            q[k+1][j]=sp.cancel(q[k][j+1]*e[k][j+1]/e[k][j])
        out.append(q[k][0]); out.append(e[k][0])
    return out
for name,c in [('ogf',P)]:
    a=qd([sp.sympify(v) for v in c],len(c)-1)
    for i,v in enumerate(a[:9]): print(i+1, sp.factor(v))
