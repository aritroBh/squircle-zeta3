import sympy as sp
x,y=sp.symbols('x y')
# clearly symmetric poly
p=x**2+y**2+x*y
print("subs test on x^2+y^2+xy:", sp.expand(p - p.subs({x:y,y:x})))  # expect 0
a1,a2,a3,a4,a5=sp.symbols('a1 a2 a3 a4 a5')
# small reversal-symmetric: a1*a5 + a2*a4 + a3^2
q=a1*a5+a2*a4+a3**2
print("subs on a1a5+a2a4+a3^2:", sp.expand(q - q.subs({a1:a5,a2:a4,a4:a2,a5:a1})))  # expect 0
# So the subs method DOES work. The earlier False on P7 must be from P7 being large -> expand timeout? re-test P7 subs
P7=sp.sympify(open('P7.txt').read())
diff=sp.expand(P7 - P7.subs({a1:a5,a2:a4,a4:a2,a5:a1}))
print("P7 subs diff == 0:", diff==0, " (len args:", len(sp.Add.make_args(diff)), ")")
