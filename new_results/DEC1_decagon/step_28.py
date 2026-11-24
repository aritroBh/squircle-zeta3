import os, sympy as sp, pickle, mpmath as mp
os.chdir('research/08_theory/DEC1_decagon')
r5 = sp.sqrt(5); phi = (1+r5)/2
t, d = sp.symbols('t d')
A0=A1=B1=1/(2*phi**2); B4=B5=1/(2*phi); kap1=phi-1; kap4=1/phi-1
h = 1/(2*phi); s36sq=(5-r5)/8; s36=sp.sqrt(s36sq); R=1/(2*phi**2*s36)
halfL2 = 1/(2*phi**6*(5-r5)); Lexpr = 2*sp.sqrt(halfL2)
al=1/(2*phi**2); Ze = 6*al*al*(h**2/Lexpr)*(3*al-2); C5 = sp.simplify(10*Ze)
pickle.dump({'A0':A0,'A1':A1,'B1':B1,'B4':B4,'B5':B5,'kap1':kap1,'kap4':kap4,'L':Lexpr,'h':h,'halfL2':halfL2,'R':R,'phi':phi,'C5':C5,'Ze':Ze}, open('geometry.pkl','wb'))
# normalized geometry
h2oL2 = (5+2*r5)/4
yypn = sp.Rational(-1,2)
yy2n = h2oL2 + (1-2*t)**2/4
L0,L1,L2,L3,L4 = sp.symbols('L0 L1 L2 L3 L4')
a  = [A0*t, A1*(1-t), B1*(1+kap1*t), B4*(1+kap4*t), B5]
ap = [sp.diff(ai,t) for ai in a]
Ls_ = [L0,L1,L2,L3,L4]
n=3
def expser(Li): return [sp.Integer(1), Li, Li**2/2, Li**3/6]
F = [sp.expand(sum(a[i]*expser(Ls_[i])[k] for i in range(5))) for k in range(n+1)]
S = [sp.expand(sum(ap[i]*expser(Ls_[i])[k] for i in range(5))) for k in range(n+1)]
print("F[0]=",sp.simplify(F[0])," S[0]=",sp.simplify(S[0]))
def pmul(P,Q):
    R=[sp.Integer(0)]*(n+1)
    for i in range(n+1):
        for j in range(n+1-i): R[i+j]=sp.expand(R[i+j]+P[i]*Q[j])
    return R
def pdiv(P,Q):
    R=[sp.Integer(0)]*(n+1)
    for i in range(n+1):
        R[i]=sp.expand((P[i]-sum(R[j]*Q[i-j] for j in range(i)))/Q[0])
    return R
def plog(P):
    Pp=[P[i]*i for i in range(1,n+1)]+[sp.Integer(0)]
    Q=pdiv(Pp,P); L=[sp.Integer(0)]
    for i in range(n): L.append(sp.expand(Q[i]/(i+1)))
    return L
def pexp(P):
    E=[sp.Integer(1)]+[sp.Integer(0)]*n
    for k in range(1,n+1): E[k]=sp.expand(sum(P[j]*j*E[k-j] for j in range(1,k+1))/k)
    return E
Phi = pdiv(S,F)
lg = plog(F)
inv = [1,-1,1,-1]
m2lg = [sp.expand(-2*x) for x in lg]
fac = pexp(pmul(m2lg,inv))
Phisq = pmul(Phi,Phi)
Q = [sp.Integer(0)]*(n+1)
Q[0] = sp.Integer(1)
for k in range(1,n+1): Q[k] = sp.expand(-2*yypn*Phi[k])
for k in range(n+1): Q[k] = sp.expand(Q[k] + yy2n*Phisq[k])
G = pmul(fac,Q)
sq=[sp.Integer(0)]*(n+1); sq[0]=sp.Integer(1)
for k in range(1,n+1): sq[k]=sp.expand((G[k]-sum(sq[j]*sq[k-j] for j in range(1,k)))/2)
E3n = sp.expand(sq[3])
print("E3n terms:", len(sp.Add.make_args(E3n)))
pickle.dump({'E3n':E3n}, open('E3n_sym.pkl','wb'))
with open('progress.log','a') as f:
    f.write(f"step: built exact normalized E3n=E3/L, {len(sp.Add.make_args(E3n))} terms, symbolic logs L0..L4. saved E3n_sym.pkl + geometry.pkl\n")
print("done")