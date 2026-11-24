# step_19.py -- assemble S_2/S_2(1) = exp(Lambda_2)*B_2, compare vs H1 Route A; build P/sqrt2 and log form.
import json
base='research/'
# Lambda_2(tau) = -sum_{k>=2} Z_k (1-2^{1-k}) tau^k/k
Lam={}
for k in range(2,KMAX+1):
    Lam[k]={(k,):-Rational(1-2**(1-k))/k}
eLam=s_exp(Lam,KMAX)
S2=s_mul(eLam,B2,KMAX)     # S_2/S_2(1)
logB2=s_log(B2,KMAX)
logS2=s_add(Lam,logB2)     # log(P/(4 sqrt2))
# --- compare with H1 routeA n=2 (eta series = S_2/S_2(1) in tau) ---
h1=json.load(open(base+'10_frontier/H1_octahedron_end/routeA_n2_K12.json'))['2']['eta']
zsyms={k:__import__('sympy').symbols(f'z{k}') for k in range(2,14)}
def to_sympy(series,k):
    out=0
    for m,c in series.get(k,{}).items():
        t=sp.Rational(c)
        for i in m: t*=zsyms[i]
        out+=t
    return sp.expand(out)
maxdiff=0
for k in range(0,13):
    h=sp.sympify(h1[k]).subs({'z%d'%i:zsyms[i] for i in range(2,14)})
    diff=sp.expand(h-to_sympy(S2,k))
    maxdiff=max(maxdiff, abs(float(diff)) if diff.is_number else (1 if diff!=0 else 0))
    if diff!=0: print("MISMATCH at k=",k, diff)
print("S_2/S_2(1) vs H1 Route A: all 13 coefficients identical:", all(sp.expand(sp.sympify(h1[k])-to_sympy(S2,k))==0 for k in range(13)))
print()
print("log(P/(4 sqrt2)) =")
show(logS2,'logS2')
print()
print("P/sqrt(2) = 4*S_2/S_2(1):")
P=s_scale(S2,4)
show(P,'P/sqrt2')
# save
def ser_to_json(S): return {str(k):{(' '.join(map(str,m))):str(c) for m,c in d.items()} for k,d in S.items()}
json.dump({'P_over_sqrt2':ser_to_json(P),'logP':ser_to_json(logS2),'B2':ser_to_json(B2),'Lambda2':{str(k):{(' '.join(map(str,m))):str(c) for m,c in d.items()} for k,d in Lam.items()}}, open('tau_series_2d.json','w'), indent=1)
print("saved tau_series_2d.json")
