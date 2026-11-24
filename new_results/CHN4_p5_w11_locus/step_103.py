import sympy as sp, pickle
a1,a2,a3,a4,a5=sp.symbols('a1 a2 a3 a4 a5')
C4=sp.sympify(open('C4.txt').read()); Q7=sp.sympify(open('Q7p.txt').read())
Q11=pickle.load(open('Q11_poly.pkl','rb'))
Q11sp=sum(sp.Rational(c.numerator,c.denominator)*sp.prod([a1,a2,a3,a4,a5][i]**e[i] for i in range(5)) for e,c in Q11.items())
# verify identity C4=Q7=Q11=0 on {a1=0,a2=a4,a5=a3}
print("Verify identities on {a1=0, a2=a4, a3=1, a5=1}:")
for tval in [1,2,3,sp.Rational(1,2),sp.Rational(7,3)]:
    sub={a1:0,a2:tval,a3:1,a4:tval,a5:1}
    print(f"  t={tval}: C4={C4.subs(sub)}, Q7'={Q7.subs(sub)}, Q11={Q11sp.subs(sub)}")
# More generally: is {a1=0, a2=a4, a5=a3} (a3 free) in V10? test with a3=2
print("\nwith a3=2, a5=2, a2=a4=t:")
for tval in [1,3]:
    sub={a1:0,a2:tval,a3:2,a4:tval,a5:2}
    print(f"  t={tval}: C4={C4.subs(sub)}, Q7'={Q7.subs(sub)}, Q11={Q11sp.subs(sub)}")
