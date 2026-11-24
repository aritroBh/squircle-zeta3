import mpmath as mp, sympy as sp, pickle
from math import factorial
from sympy.functions.combinatorial.numbers import stirling
mp.mp.dps=90
d_,m_=sp.symbols('d m'); I=sp.I
A={dl:sp.lambdify((d_,m_),sp.sympify(pickle.load(open(f'../../a{dl}_closed.pkl','rb')),locals={'d':d_,'m':m_,'I':I}),'mpmath') for dl in (9,10)}
a=1+mp.sqrt(2); b=mp.sqrt(2)-1
def Afun(dd,p):
    q=dd-p; D=2*dd
    def S(t):
        tot=mp.mpf(0)
        for j in range(1,dd+1):
            Cj=mp.binomial(mp.mpf(1)/2,j)*mp.factorial(j)*2**dd
            for s in range(max(0,j-q),min(j,p)+1):
                tot+=Cj*int(stirling(p,s,kind=2))*int(stirling(q,j-s,kind=2))*t**(2*(j-s))*(1-t)**(2*s)*(t**2+(1-t)**2)**(dd-j)
        return tot
    def f(v):
        t=(v*v+2*v-1)/(4*v)
        return v**(2*dd-1-D)*(4*v)**D*S(t)*(v+1j*0+mp.mpc(0,1))**(2-2*dd)/((v*v-a*a)*(v*v-b*b))
    return f
for dd,p in ((11,6),(12,7),(13,7),(14,8),(16,9),(17,9)):
    f=Afun(dd,p)
    Np=128; r0=mp.mpf(3)/10; c0=mp.mpc(0,1)
    vals=[f(c0+r0*mp.expjpi(mp.mpf(2*k)/Np)) for k in range(Np)]
    tc=[sum(vals[k]*mp.expjpi(-mp.mpf(2*k*n)/Np) for k in range(Np))/Np/r0**n for n in range(11)]
    row=[]
    for dl in (9,10):
        r=tc[dl]/tc[0]
        e=A[dl](mp.mpf(dd),mp.mpf(p))
        row.append(mp.nstr(abs(r-e)/abs(e),5))
    print(dd,p,'rel err a_9,a_10 vs direct Taylor of A(v):',row,flush=True)
