import os, time
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
t0=time.time()
import sympy as sp
x,y,z,w=sp.symbols('x y z w')
C3 = (4*w**3 + 3*w**2*x + 6*w**2*y + 6*w**2*z + 3*w*x**2 - 6*w*x*y - 6*w*x*z - 6*w*y**2 + 6*w*y*z - 6*w*z**2 + 4*x**3 + 6*x**2*y + 6*x**2*z - 6*x*y**2 + 6*x*y*z - 6*x*z**2 - 4*y**3 - 3*y**2*z - 3*y*z**2 - 4*z**3)
s,t=sp.symbols('s t')
tests = {
 'Pi1 {a1=a2,a3=a4}': C3.subs({x:s,y:s,z:t,w:t}),
 'Pi2 {a1=a3,a2=a4}': C3.subs({x:s,y:t,z:s,w:t}),
 'palindromic (s,t,t,s)': C3.subs({x:s,y:t,z:t,w:s}),
 'plane a1+a4=a2+a3 (w=y+z-x)': C3.subs({w:y+z-x}),
 'locus a1=a2=a3 (s,s,s,t)': C3.subs({x:s,y:s,z:s,w:t}),
}
for name,val in tests.items():
    print(('VANISHES ' if sp.expand(val)==0 else 'NONZERO  ')+name, '' if sp.expand(val)==0 else sp.factor(val))
print('C3(1,1,1,1) =', C3.subs({x:1,y:1,z:1,w:1}))
grad=[sp.diff(C3,v).subs({x:1,y:1,z:1,w:1}) for v in (x,y,z,w)]
print('grad C3 at (1,1,1,1):', grad, ' prop. to (1,-1,-1,1)?', [sp.Rational(g)/grad[0] for g in grad])
# smoothness / absolute irreducibility of C3: groebner of (C3, partials)
t=time.time()
I=sp.groebner([C3,sp.diff(C3,x),sp.diff(C3,y),sp.diff(C3,z),sp.diff(C3,w)], x,y,z,w)
print('groebner basis of (C3,grad):', I.polys, ' (%.0fs)'%(time.time()-t))
print('t=%.1fs'%(time.time()-t0))
