import json, mpmath as mp
mp.mp.dps=40
z=mp.exp(2j*mp.pi/8)
row=json.load(open('research/08_theory/A3f_level8_relations/weight2_new_relation.json'))['row']
f=lambda j,t: z**j/(1-z**j*t)
H1=lambda j,t: -mp.log(1-z**j*t)
def H2(a,b):  # first index outermost: int_0^1 f_a(t) H(b;t) dt
    return mp.quad(lambda t: f(a,t)*H1(b,t),[0,0.5,1])
def H2r(a,b): # reversed: first index innermost
    return H2(b,a)
for name,H in [('outer-first',H2),('inner-first',H2r)]:
    try:
        s=mp.fsum(c*H(w[0],w[1]) for w,c in row); print(name,'residual',mp.nstr(abs(s),5))
    except Exception as e: print(name,'error',e)
