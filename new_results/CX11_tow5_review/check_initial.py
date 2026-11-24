"""Independent Bell recurrence, tau definitions, and freshly generated TOWN1 towers."""
import pickle
from math import comb
from fractions import Fraction as Fr
import sympy as sp
import tow5lib as T5, tow2lib as T2
import trunc_engine as te
from fastred import QSI
from driver import q4
D,M=T5.d,T5.m
for n in (9,10):
 for fn,ref in [(T2.tau1,T5.tau1f),(T2.tau2,T5.tau2f)]:
  assert sp.simplify(sp.expand_complex(fn(n))-ref(n))==0
 assert sp.conjugate(T5.tau1f(n))==(-1)**n*T5.tau2f(n)
 print('tau definitions + conjugation PASS',n,T5.tau1f(n),T5.tau2f(n),flush=True)
a=[sp.Integer(1)]
for n in range(1,11):
 a.append(sp.cancel(sum(r*T5._L[r]*a[n-r] for r in range(1,n+1))/n))
 assert sp.cancel(a[n]-T5.arho(n))==0
 print('Independent Bell recurrence vs partition assembly PASS',n,flush=True)
for dl in (9,10):
 raw=pickle.load(open(f'a{dl}_closed.pkl','rb'))
 print('closed pickle type',dl,type(raw).__name__,flush=True)
# Separate coefficient recurrence for log S, up to fifth cumulant.
r=[sp.Integer(1)]+[T2.R[j] for j in range(1,5)]+[T5.R5]
kappa={}
for n in range(1,6):
 kappa[n]=sp.cancel(r[n]-sum(comb(n-1,j-1)*kappa[j]*r[n-j] for j in range(1,n)))
assert sp.cancel(kappa[5]-T5.sl5)==0
print('Independent fifth cumulant recurrence PASS',flush=True)
# Use the freshly rebuilt cache, which is separately compared against the submission.
S2=te.stirling_table(18);fresh={};counts=0
for dd in (7,8,11,13):
 LS=te.level_series(dd,10,S2)
 for pp in sorted({dd,dd-1,(dd+1)//2}):
  vals=te.rho_key(LS,pp,dd-pp,S2)
  fresh[(dd,pp)]=[q4(v) for v in vals]
  for dl in range(11):
   v=sp.expand(a[dl].subs({D:dd,M:pp}))
   q=QSI(Fr(sp.re(v)),0,Fr(sp.im(v)),0)
   assert vals[dl]==vals[0]*q,(dd,pp,dl)
   counts+=1
 print('Fresh initial series all depths PASS d=',dd,flush=True)
print('Fresh initial coefficient equalities:',counts,flush=True)
# At two levels, propagate the independently generated initial Laurent series.
for k in (14,16):
 rho=te.rho_level(k,10,S2);tower=te.build_tower_level(k,10,rho)
 for key in [(7,0),(4,3),(5,3),(4,4),(7,4)]:
  pp,qq=key;dd=pp+qq
  base=te.rho_key(te.level_series(dd,10,S2),pp,qq,S2)[0]
  for dl in (9,10):
   if pp==qq:
    import prove_E3b910 as pe
    coeff=[pe.ehat(dl,j).subs(pe.p,pp) for j in range(dl//2+1)]
    z=k-dd-1
   else:
    coeff=[T5.forced(dl,j,dd,pp) for j in range(dl//2+1)];z=k-dd
   v=sp.cancel(sum(comb(z,j)*coeff[j] for j in range(dl//2+1)))
   pred=base*QSI(Fr(sp.re(v)),0,Fr(sp.im(v)),0)
   assert tower[key][dl]==pred,(k,key,dl)
 print('Fresh propagated tower: 10 depth-9/10 cells PASS k=',k,flush=True)
pickle.dump(fresh,open('fresh_initial.pkl','wb'))
