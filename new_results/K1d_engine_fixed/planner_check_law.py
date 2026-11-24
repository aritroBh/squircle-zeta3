import sympy as sp, json, mpmath as mp
R=sp.Rational
def Phi(A0,A1,B,r,L2,q0,V02):
    k=r-1; c=(q0**2-L2*V02)/L2**R(5,2)
    PA=6*B*c*(3*A0*k*q0-L2*(2*A0+r*A1))
    PL=18*B*c*(k*(-A1**2+A1*B*(2*k-1)+2*B**2*k)*q0+L2*(A0*B+r*A1**2))
    return PA,PL
cont={2:(R(-5,8),R(-3,4)), R(1,2):(1,R(7,8))}
cases=[('b0 e1',R(1,2),R(3,4),R(1,4),2,5,1,1,R(159,500)*sp.sqrt(5)),
       ('b0 e2',R(3,4),R(1,2),R(1,2),R(1,2),5,-6,8,R(159,500)*sp.sqrt(5)),
       ('b1 e1',R(1,2),R(3,4),R(1,4),2,13,3,1,R(417,8788)*sp.sqrt(13)),
       ('b1 e2',R(3,4),R(1,2),R(1,2),R(1,2),10,-14,20,R(21,250)*sp.sqrt(10))]
for name,A0,A1,B,r,L2,q0,V02,claim in cases:
    PA,PL=Phi(A0,A1,B,sp.nsimplify(r),L2,q0,V02); a,l=cont[sp.nsimplify(r)]
    Ne=sp.nsimplify(sp.simplify(PA*a+PL*l)); print(name,'N_e =',Ne,'  claim',claim,'  match',sp.simplify(Ne-claim)==0)
# independent check of the b=1 PSLQ relations on the verifier's 110-digit data
mp.mp.dps=110
d=json.load(open('research/10_frontier/E1_edge_law/verify/out/t3_b1_110.json'))
z3=mp.zeta(3);L=mp.log(2);pi=mp.pi
rels={1:[8788,-615,-297,4748,208,16848,-40560,0],2:[2000,-267,-132,3336,80,-2160,-7200,0],3:[800,15,-4,-6400,4,19200,-9600,0]}
for a in d['arcs']:
    e=int(a['edge']); X=mp.mpf(a['P3_arc'])/mp.mpf(a['P0'])
    v=[X,z3,pi**2*L,L**3,pi**2,L**2,L,1]
    print('b1 edge',e,'residual',mp.nstr(mp.fsum(c*x for c,x in zip(rels[e],v)),3))
