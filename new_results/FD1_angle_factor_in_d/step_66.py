import sympy as sp, math
y,x=sp.symbols('y x')
P={0:sp.Integer(1),1:sp.Integer(1),2:1-x,3:2-x,4:x**2-6*x+5,5:x**2-13*x+16,6:-x**3+29*x**2-89*x+61,
   7:-x**3+60*x**2-297*x+272,8:x**4-124*x**3+1134*x**2-2396*x+1385,
   9:x**4-251*x**3+3651*x**2-10841*x+7936,10:-x**5+507*x**4-12212*x**3+62372*x**2-101187*x+50521}
N=9
G=sp.expand(sum(P[d]*y**d/math.factorial(d) for d in range(11)))
Gm=sp.expand(G.subs(y,-y))
def ser(expr, n=N):
    return sp.series(expr, y, 0, n).removeO().expand()
prod=ser(G*Gm)
print("G(y)G(-y) coeffs (y^d/d!):")
for d in range(0,N):
    c=sp.expand(prod).coeff(y,d)*math.factorial(d)
    print(d, sp.factor(c))
