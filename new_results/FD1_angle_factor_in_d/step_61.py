import sympy as sp
from fractions import Fraction as Fr
x=sp.symbols('x')
P={2:1-x,3:2-x,4:x**2-6*x+5,5:x**2-13*x+16,6:-x**3+29*x**2-89*x+61,
   7:-x**3+60*x**2-297*x+272,8:x**4-124*x**3+1134*x**2-2396*x+1385,
   9:x**4-251*x**3+3651*x**2-10841*x+7936,10:-x**5+507*x**4-12212*x**3+62372*x**2-101187*x+50521}
# include P_0=1, P_1=1
P[0]=sp.Integer(1); P[1]=sp.Integer(1)
MAXD=5
def vec(p):
    return [sp.expand(p).coeff(x,i) for i in range(MAXD+1)]
# Try ODE: P_{d+r} = sum_{i=0}^{r-1} a_i(x) P_{d+i}, a_i polynomial deg<=Dg
for r in [1,2,3]:
    for Dg in [0,1,2]:
        # unknowns: coeffs of a_i
        nunk = r*(Dg+1)
        unk = sp.symbols('u0:'+str(nunk))
        eqs=[]
        for d in range(0, 11-r):
            lhs = vec(P[d+r])
            rhs = [0]*(MAXD+1)
            for i in range(r):
                ai = sum(unk[i*(Dg+1)+k]*x**k for k in range(Dg+1))
                term = sp.expand(ai*P[d+i])
                rhs = [rhs[j]+sp.expand(term).coeff(x,j) for j in range(MAXD+1)]
            for j in range(MAXD+1):
                eqs.append(sp.Eq(lhs[j], rhs[j]))
        sol = sp.solve(eqs, unk, dict=True)
        if sol:
            print(f"r={r}, Dg={Dg}: SOLUTION FOUND: {sol[0]}")
        else:
            print(f"r={r}, Dg={Dg}: none")
