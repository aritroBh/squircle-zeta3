import os, mpmath as mp, itertools, pickle, math
os.chdir('research/08_theory/DEC1_decagon')
mp.mp.dps = 100
s5 = mp.sqrt(5); phi = (1+s5)/2
A0=1/(2*phi**2); A1=A0; B1=A0; B4=1/(2*phi); B5=B4
kap1=phi-1; kap4=1/phi-1
h2oL2=(5+2*s5)/4
KEYS = [k for k in itertools.product(range(4),repeat=5) if sum(k)<=3]
def eadd(P,Q):
    R={}
    for k in set(P)|set(Q):
        aa=P.get(k,[mp.mpf(0)]*4); b=Q.get(k,[mp.mpf(0)]*4)
        R[k]=[aa[i]+b[i] for i in range(4)]
    return R
def ekmul(P,Q):
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
def escale(P,s): return {k:[v[i]*s for i in range(4)] for k,v in P.items()}
def ediv_scalar(P,s): return {k:[v[i]/s for i in range(4)] for k,v in P.items()}
ONE={(0,0,0,0,0):[mp.mpf(1),0,0,0]}
def pmul(P,Q):
    R=[{} for _ in range(4)]
    for i in range(4):
        for j in range(4-i): R[i+j]=eadd(R[i+j],ekmul(P[i],Q[j]))
    return R
def pdiv(P,Q):
    R=[{} for _ in range(4)]
    q0=Q[0][(0,0,0,0,0)][0]
    for i in range(4):
        acc=P[i]
        for j in range(i): acc=eadd(acc,escale(ekmul(R[j],Q[i-j]),mp.mpf(-1)))
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
apol = [[mp.mpf(0),A0],[A1,-A1],[B1,B1*kap1],[B4,B4*kap4],[B5,mp.mpf(0)]]
appol = [[A0,mp.mpf(0)],[-A1,mp.mpf(0)],[B1*kap1,mp.mpf(0)],[B4*kap4,mp.mpf(0)],[mp.mpf(0),mp.mpf(0)]]
F=[{} for _ in range(4)]; S=[{} for _ in range(4)]
for k in range(4):
    fac=mp.mpf(1)/math.factorial(k)
    for i in range(5):
        for e in range(5):
            key=tuple(k if j==e else 0 for j in range(5))
            F[k]=eadd(F[k],{key:[apol[i][0]*fac, apol[i][1]*fac,0,0]})
            S[k]=eadd(S[k],{key:[appol[i][0]*fac, appol[i][1]*fac,0,0]})
print("F0:",F[0]); print("S0:",S[0])
Phi=pdiv(S,F)
lg=plog(F)
inv=[mp.mpf(1),mp.mpf(-1),mp.mpf(1),mp.mpf(-1)]
m2lg=[escale(x,mp.mpf(-2)) for x in lg]
facE=pexp(pmul(m2lg,inv))
Phisq=pmul(Phi,Phi)
Q=[{} for _ in range(4)]
Q[0]=ONE
for k in range(1,4): Q[k]=escale(Phi[k], mp.mpf(1))  # -2*yypn = 1
yy2poly=[h2oL2+mp.mpf('0.25'), mp.mpf(-1), mp.mpf(1), mp.mpf(0)]
for k in range(4):
    Q[k]=eadd(Q[k], ekmul(Phisq[k], {(0,0,0,0,0):yy2poly}))
G=pmul(facE,Q)
sq=[{} for _ in range(4)]; sq[0]=ONE
for k in range(1,4):
    acc=G[k]
    for j in range(1,k): acc=eadd(acc,escale(ekmul(sq[j],sq[k-j]),mp.mpf(-1)))
    sq[k]=escale(acc,mp.mpf('0.5'))
E3n=sq[3]
nz=[k for k in E3n if any(abs(x)>mp.mpf('1e-80') for x in E3n[k])]
print("E3n monomials present:", len(nz))
pickle.dump(E3n, open('E3n_num100.pkl','wb'))
with open('progress.log','a') as f:
    f.write(f"step: numeric (100dps) series engine -> E3n_num100.pkl; {len(nz)} nonzero keys\n")
print("saved")