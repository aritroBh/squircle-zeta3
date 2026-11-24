# step_29.py -- direct perimeter quadrature (F1 support-function form + arclength form) at p=1.02,1.05, 60 dps;
# truncation table of the tau-series (my 64 numeric coeffs); agreement digits.
import mpmath as mp, json
mp.mp.dps=70
S=[mp.mpf(x) for x in json.load(open('hi_tau_numeric_n2.json'))]
def P_quad(p):
    q=p/(p-1)
    f=lambda u: (1+u**q)**(1/q)*(1+u**2)**(-mp.mpf(3)/2)
    return 8*mp.quad(f,[0,1])
def P_arc(p):
    # x=cos^{2/p} t, y=sin^{2/p} t ; dx,dy; integrand sqrt(x'^2+y'^2)
    e=2/p
    f=lambda t: mp.sqrt((e*mp.cos(t)**(e-1)*mp.sin(t))**2+(e*mp.sin(t)**(e-1)*mp.cos(t))**2)
    return 4*mp.quad(f,[0,mp.pi/2])
res={}
for p in [mp.mpf('1.02'),mp.mpf('1.05')]:
    P1=P_quad(p); P2=P_arc(p)
    print(f"p={p}: support-function quad P={mp.nstr(P1,50)}")
    print(f"          arclength quad   P={mp.nstr(P2,50)}  |diff|={mp.nstr(abs(P1-P2),3)}")
    res[str(p)]=P1
tau={ '1.02': 1-1/mp.mpf('1.02'), '1.05': 1-1/mp.mpf('1.05')}
print()
for ps in ['1.02','1.05']:
    t=tau[ps]; P=res[ps]
    print(f"--- p={ps}, tau={mp.nstr(t,12)} ---")
    print(" N    trunc. abs error      est. next-term")
    terms=[S[k]*t**k for k in range(65)]
    for Nn in [4,8,12,16,20,24,32,40,48,56,64]:
        T=4*mp.sqrt(2)*sum(terms[:Nn+1])
        err=abs(T-P)
        print(f"{Nn:3d}  {mp.nstr(err,3):>12}   {mp.nstr(abs(terms[Nn+1]) if Nn+1<=64 else 0,3):>12}")
    # least term
    mags=[(abs(terms[k]),k) for k in range(1,65)]
    mags.sort()
    print(" least |term| among k<=64:", mp.nstr(mags[0][0],3), "at k=",mags[0][1])
json.dump({k:mp.nstr(v,60) for k,v in res.items()}, open('quad_ref.json','w'))
