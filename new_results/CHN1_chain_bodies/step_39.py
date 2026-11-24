import os, time, json
WORK='<project>/research/08_theory/CHN1_chain_bodies'
os.chdir(WORK)
print('blind points so far:', sum(1 for _ in open('blind_points.jsonl')) if os.path.exists('blind_points.jsonl') else 0)
t0=time.time()
import sympy as sp
x,y,z,w,s,t=sp.symbols('x y z w s t')
Q5 = (64*w**5 + 80*w**4*x + 160*w**4*y + 160*w**4*z + 61*w**3*x**2 + 64*w**3*x*y + 64*w**3*x*z + 64*w**3*y**2 + 320*w**3*y*z + 64*w**3*z**2 + 61*w**2*x**3 + 90*w**2*x**2*y + 90*w**2*x**2*z - 96*w**2*x*y**2 + 93*w**2*x*y*z - 96*w**2*x*z**2 - 64*w**2*y**3 + 96*w**2*y**2*z + 96*w**2*y*z**2 - 64*w**2*z**3 + 80*w*x**4 + 64*w*x**3*y + 64*w*x**3*z - 96*w*x**2*y**2 + 93*w*x**2*y*z - 96*w*x**2*z**2 - 320*w*x*y**3 - 93*w*x*y**2*z - 93*w*x*y*z**2 - 320*w*x*z**3 - 160*w*y**4 - 64*w*y**3*z - 90*w*y**2*z**2 - 64*w*y*z**3 - 160*w*z**4 + 64*x**5 + 160*x**4*y + 160*x**4*z + 64*x**3*y**2 + 320*x**3*y*z + 64*x**3*z**2 - 64*x**2*y**3 + 96*x**2*y**2*z + 96*x**2*y**2*z + 96*x**2*y*z**2 - 64*x**2*z**3 - 160*x*y**4 - 64*x*y**3*z - 90*x*y**2*z**2 - 64*x*y*z**3 - 160*x*z**4 - 64*y**5 - 80*y**4*z - 61*y**3*z**2 - 61*y**2*z**3 - 80*y*z**4 - 64*z**5)
# NOTE: duplicated 96*x**2*y**2*z term by accident? check factor below vs original
print('Q5(1,1,1,1) =', Q5.subs({x:1,y:1,z:1,w:1}))
for name,val in {'Pi1':Q5.subs({x:s,y:s,z:t,w:t}),'Pi2':Q5.subs({x:s,y:t,z:s,w:t}),
                 'a1+a4=a2+a3':Q5.subs({w:y+z-x}),'pal (s,t,t,s)':Q5.subs({x:s,y:t,z:t,w:s})}.items():
    val=sp.expand(val)
    print(('VANISHES ' if val==0 else 'NONZERO  ')+name, '' if val==0 else sp.factor(val))
# middle factor of z335 form
mid = w**2*x + w*x**2 + 2*w*x*y + 2*w*x*z + 2*w*y*z + 2*x*y*z + y**2*z + y*z**2
print('mid factor:', sp.factor(mid))
print('t=%.1fs'%(time.time()-t0))
