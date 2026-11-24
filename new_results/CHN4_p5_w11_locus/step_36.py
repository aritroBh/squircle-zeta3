import pickle
from fractions import Fraction as F
import sympy as sp
Q11=pickle.load(open('Q11_poly.pkl','rb'))
# 1. reversal symmetry
revok=all(Q11.get(e,F(0))==Q11.get(e[::-1],F(0)) for e in set(list(Q11)+[e[::-1] for e in Q11]))
print("Q11 reversal-symmetric:", revok)
# homogeneity
degs={sum(e) for e in Q11}
print("homogeneous degrees present:", degs)

# 2. build sympy expr
a1,a2,a3,a4,a5=sp.symbols('a1 a2 a3 a4 a5')
A=[a1,a2,a3,a4,a5]
Q11sp=sum(sp.Rational(c.numerator,c.denominator)*sp.prod(A[i]**e[i] for i in range(5)) for e,c in Q11.items())
print("Q11 sympy terms:", len(sp.Add.make_args(Q11sp)))

# 3. evaluate at planner V10 points (numerically) and compare f11/L
import mpmath as mp
mp.mp.dps=50
f=sp.lambdify(A, Q11sp, 'mpmath')
pts=[(mp.mpf('1.29366878462024308244448284389'), mp.mpf('0.070519999391807054958080626713')),
     (mp.mpf('5.19669194677918082507395569035'), mp.mpf('6.53040365342565411592415935716'))]
plan_f11=[mp.mpf('-4.2895814'), mp.mpf('-10747354.0')]
for (a4v,a5v),pf in zip(pts,plan_f11):
    qv=f(1,1,1,a4v,a5v)
    L=1*1*a4v*(1+1+1+a4v+a5v)
    print(f"  V10 pt a4={mp.nstr(a4v,8)}: Q11={mp.nstr(qv,10)}, Q11*L={mp.nstr(qv*L,10)} vs planner f11~{mp.nstr(pf,10)}")

# 4. factor over Q
print("\nfactoring Q11 over Q...")
t0=__import__('time').time()
fac=sp.factor(Q11sp)
print("factor time:", __import__('time').time()-t0)
print("factor structure:", fac)
