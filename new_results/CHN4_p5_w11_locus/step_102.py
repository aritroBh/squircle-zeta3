import sympy as sp, mpmath as mp, pickle, time
mp.mp.dps=60
a1,a2,a3,a4,a5=sp.symbols('a1 a2 a3 a4 a5')
C4=sp.sympify(open('C4.txt').read()); Q7=sp.sympify(open('Q7p.txt').read())
Q11=pickle.load(open('Q11_poly.pkl','rb'))
Q11sp=sum(sp.Rational(c.numerator,c.denominator)*sp.prod([a1,a2,a3,a4,a5][i]**e[i] for i in range(5)) for e,c in Q11.items())
t=sp.symbols('t')
# impose a1=0,a3=1,a5=1,a2=a4=t
sub={a1:0,a2:t,a3:1,a4:t,a5:1}
c4t=sp.expand(C4.subs(sub)); q7t=sp.expand(Q7.subs(sub)); q11t=sp.expand(Q11sp.subs(sub))
print("C4(0,t,1,t,1) =", sp.factor(c4t))
print("Q7(0,t,1,t,1) =", sp.factor(q7t))
print("Q11(0,t,1,t,1)=", sp.factor(q11t))
# common roots
import numpy as np
for nm,poly in [('C4',c4t),('Q7',q7t),('Q11',q11t)]:
    rr=sp.nroots(poly, n=30)
    pos=[sp.re(r) for r in rr if abs(sp.im(r))<sp.Float('1e-20') and sp.re(r)>0]
    print(f"{nm} positive roots:", [mp.nstr(mp.mpf(str(x)),20) for x in pos])
def log_progress(line):
    with open('progress.log','a') as f: f.write(time.strftime('%H:%M:%S ')+line+'\n')
log_progress("Final W12 batch DONE: 8 Gamma pts (a1=5e-5..2.0), all controls vanish, all W12 irreducibles NONZERO. Normalized ratio max|irred|/max|c12| always positive: 6.1e-8 (a1=5e-5) up to 5.6e-3 (a1=1). Ratio ~ LINEAR in a1 as a1->0 => W12 irreducibles -> 0 at boundary pt (0,2.667,1,2.667,1) [a1=0, NOT open orthant]. Investigating boundary point directly.")
