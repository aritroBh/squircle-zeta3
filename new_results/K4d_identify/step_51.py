import mpmath as mp, json
from mpmath import mpf
mp.mp.dps=80
WORK='research/08_theory/K4d_identify'
def load(tag): return [mpf(s) for s in json.load(open(WORK+'/modelA_%s.json'%tag))['Sdk']]
S64,S96,S128,S192 = load('h64'),load('h96'),load('h128'),load('h192')
def digs(a,b): 
    d=abs(a-b); return float(-mp.log10(d)) if d>0 else 200.0
print("Model A digit agreements:")
for k in range(5):
    print(" d^%d: h64/96=%.1f  h96/128=%.1f  h128/192=%.1f"%(k,digs(S64[k],S96[k]),digs(S96[k],S128[k]),digs(S128[k],S192[k])))
# Richardson: err ~ C*h^p? tanh-sinh err ~ exp(-c/sqrt(h)) i.e. log err ~ -c/sqrt(h)
# use three finest: h=1/96,1/128,1/192 (h96,h128,h192)
import math
print("\nRichardson/shanks extrapolation (using h96,h128,h192):")
best={}
for k in range(5):
    # Aitken / assume err_k = A*r^(1/h)? Use ratio of diffs
    d1=S128[k]-S96[k]; d2=S192[k]-S128[k]
    if d1!=0 and d2!=0 and d1!=d2:
        # assume S(h)=S*+C*q^(1/h): 1/h = 96,128,192. ratio (S128-S96)/(S192-S128)=q^(96-128)... 
        # S96=S*+C q^96, S128=S*+C q^128, S192=S*+C q^192
        # d1=C(q^128-q^96)=C q^96(q^32-1); d2=C q^128(q^64-1)
        # d2/d1 = q^32 (q^64-1)/(q^32-1) = q^32(q^32+1)
        r=d2/d1
        # solve q^32(q^32+1)=r  => y^2+y-r=0, y=q^32
        y=(-1+mp.sqrt(1+4*r))/2
        if y>0 and y<1:
            q=y**(mpf(1)/32)
            Sstar=S192[k] + S192[k]*0  # placeholder
            # S* = S192 + C q^192 ; C = d2/(q^128(q^64-1))
            C=d2/(q**128*(q**64-1))
            Sstar=S192[k]-C*q**192
            err=abs(C*q**192)
            best[k]=Sstar
            print(" d^%d: q=%.4f extrap err~%.1e  S*=%s"%(k,float(q),float(err),mp.nstr(Sstar,20)))
        else:
            best[k]=S192[k]; print(" d^%d: y=%.3g not in (0,1), use h192"%(k,float(y)))
print("\nh192 values:")
for k in range(5): print("  d^%d = %s"%(k,mp.nstr(S192[k],60)))
