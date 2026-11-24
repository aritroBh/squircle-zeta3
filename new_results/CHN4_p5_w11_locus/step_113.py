import sympy as sp
x,y=sp.symbols('x y')
print("sequential-check x^2+y^2+xy via simultaneous=True:", sp.expand((x**2+y**2+x*y) - (x**2+y**2+x*y).subs({x:y,y:x}, simultaneous=True)))
a1,a2,a3,a4,a5=sp.symbols('a1 a2 a3 a4 a5')
P7=sp.sympify(open('P7.txt').read())
d=sp.expand(P7 - P7.subs({a1:a5,a2:a4,a4:a2,a5:a1}, simultaneous=True))
print("P7 reversal-symmetric (simultaneous=True):", d==0)
# Q7', C4 reversal check too (should be symmetric per CHN3)
for nm in ['C4.txt','Q7p.txt']:
    e=sp.sympify(open(nm).read())
    dd=sp.expand(e - e.subs({a1:a5,a2:a4,a4:a2,a5:a1}, simultaneous=True))
    print(f"{nm} reversal-symmetric:", dd==0)
