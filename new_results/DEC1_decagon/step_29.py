import os, mpmath as mp, itertools, pickle
os.chdir('research/08_theory/DEC1_decagon')
mp.mp.dps = 100
s5 = mp.sqrt(5); phi = (1+s5)/2
A0=1/(2*phi**2); A1=A0; B1=A0; B4=1/(2*phi); B5=B4
kap1=phi-1; kap4=1/phi-1
h2oL2=(5+2*s5)/4
n=3
# monomial keys: weight<=3 in 5 vars
KEYS = [k for k in itertools.product(range(4),repeat=5) if sum(k)<=3]
KIDX = {k:i for i,k in enumerate(KEYS)}
NK = len(KEYS)
print("log monomials:", NK)
# element = dict key -> [c0,c1,c2,c3] mpf
def eadd(P,Q):
    R={}
    for k in set(P)|set(Q):
        a=P.get(k,[mp.mpf(0)]*4); b=Q.get(k,[mp.mpf(0)]*4)
        R[k]=[a[i]+b[i] for i in range(4)]
    return R
def ekmul(P,Q):  # multiply, truncate weight 3
    R={}
    for k1,v1 in P.items():
        for k2,v2 in Q.items():
            k=tuple(k1[i]+k2[i] for i in range(5))
            if sum(k)>3: continue
            c=[mp.mpf(0)]*4
            for i in range(4):
                for j in range(4-i): c[i+j]+=v1[i]*v2[j]
            if k in R: R[k]=[R[k][i]+c[i] for i in range(4)]
            else: R[k]=c
    return R
def escale(P,s):
    return {k:[v[i]*s for i in range(4)] for k,v in P.items()}
def ediv_scalar(P,s):
    return {k:[v[i]/s for i in range(4)] for k,v in P.items()}
ZERO={}; ONE={(0,0,0,0,0):[mp.mpf(1),0,0,0]}
# series in d: list of 4 elements
def pmul(P,Q):
    R=[{} for _ in range(4)]
    for i in range(4):
        for j in range(4-i):
            R[i+j]=eadd(R[i+j],ekmul(P[i],Q[j]))
    return R
def pdiv(P,Q):
    R=[{} for _ in range(4)]
    q0=Q[0][(0,0,0,0,0)][0]
    for i in range(4):
        acc=P[i]
        for j in range(i):
            acc=eadd(acc,escale(ekmul(R[j],Q[i-j]),mp.mpf(-1)))
        R[i]=ediv_scalar(acc,q0)
    return R
def plog(P):
    Pp=[escale(P[i],mp.mpf(i)) for i in range(1,4)]+[{}]
    Q=pdiv(Pp,P)
    L=[{}]
    for i in range(3): L.append(escale(Q[i],mp.mpf(1)/(i+1)))
    return L
def pexp(P):
    E=[ONE]+[{}]*3
    for k in range(1,4):
        acc={}
        for j in range(1,k+1): acc=eadd(acc,escale(ekmul(P[j],E[k-j]),mp.mpf(j)))
        E[k]=escale(acc,mp.mpf(1)/k)
    return E
a  = [lambda tt: A0*tt, lambda tt: A1*(1-tt), lambda tt: B1*(1+kap1*tt), lambda tt: B4*(1+kap4*tt), lambda tt: B5]
apol = [[mp.mpf(0),A0],[A1,-A1],[B1,B1*kap1],[B4,B4*kap4],[B5,mp.mpf(0)]]
appol = [[A0,0],[-A1,0],[B1*kap1,0],[B4*kap4,0],[mp.mpf(0),0]]
F=[{} for _ in range(4)]; S=[{} for _ in range(4)]
import math
for k in range(4):
    for i in range(5):
        for e in range(5):
            key=tuple(k if j==e else 0 for j in range(5))
            c=escale({key:[apol[i][0],apol[i][1] if k>=0 else 0]}, mp.mpf(1)/math.factorial(k))
            # note t-coeff: a_i(t)=apol[i][0]+apol[i][1] t
            c={key:[apol[i][0]/math.factorial(k), apol[i][1]/math.factorial(k),0,0]}
            F[k]=eadd(F[k],c)
            c2={key:[appol[i][0]/math.factorial(k), appol[i][1]/math.factorial(k),0,0]}
            S[k]=eadd(S[k],c2)
print("F0:",F[0], "S0:",S[0])
Phi=pdiv(S,F)
lg=plog(F)
inv=[mp.mpf(1),mp.mpf(-1),mp.mpf(1),mp.mpf(-1)]
m2lg=[escale(x,mp.mpf(-2)) for x in lg]
fac=pexp(pmul(m2lg,inv))
Phisq=pmul(Phi,Phi)
Q=[{} for _ in range(4)]
Q[0]=ONE
for k in range(1,4): Q[k]=escale(Phi[k], -2*(-mp.mpf('0.5')))  # -2*yypn*Phi[k], yypn=-1/2
for k in range(4):
    # + yy2n * Phisq[k]; yy2n = h2oL2 + (1-2t)^2/4 -> t-poly [h2oL2+1/4, -1, 1]
    yy2poly=[h2oL2+mp.mpf('0.25'), mp.mpf(-1), mp.mpf(1), mp.mpf(0)]
    Q[k]=eadd(Q[k], ekmul(Phisq[k], {(0,0,0,0,0):yy2poly}))
G=pmul(fac,Q)
sq=[{} for _ in range(4)]; sq[0]=ONE
for k in range(1,4):
    acc=G[k]
    for j in range(1,k):
        acc=eadd(acc,escale(ekmul(sq[j],sq[k-j]),mp.mpf(-1)))
    sq[k]=escale(acc,mp.mpf('0.5'))
E3n=sq[3]
print("E3n monomials present:", len([k for k in E3n if any(abs(x)>mp.mpf('1e-90') for x in E3n[k])]))
pickle.dump(E3n, open('E3n_num100.pkl','wb'))
with open('progress.log','a') as f:
    f.write(f"step: numeric (100dps) series engine -> E3n_num100.pkl; {len(E3n)} keys\n")
print("saved")