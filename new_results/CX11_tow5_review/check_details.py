"""Finite master identities, exact prefactors, printed formulas, and seed nonuniqueness."""
from pathlib import Path
import pickle,re
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication_application, convert_xor
import tow5lib as T, tow2lib as T2, tow3lib as T3
import prove_E3b910 as pe
D,M=T.d,T.m;I=sp.I
local={'d':D,'m':M,'i':I,'Pi_4':T.Pi(4),'Pi_5':T.Pi(5)}
trans=standard_transformations+(implicit_multiplication_application,convert_xor)
def parse(x): return parse_expr(x.replace('\n',' '),local_dict=local,transformations=trans)
proof=Path('../TOW5_tower_depth10/proof_depth10.md').read_text()
for dl,ss,den in [(9,2,240),(9,3,12),(10,3,96)]:
 match=re.search(rf'  F\({dl},{ss}\) = \((.*?)\) / \({den} Pi_{dl//2}\)',proof,re.S);assert match
 assert sp.cancel(parse(match[1])/(den*T.Pi(dl//2))-T.forced(dl,ss))==0
 print('Printed F formula PASS',dl,ss,flush=True)
for dl,ss in [(9,4),(10,4),(10,5)]:
 match=re.search(rf'  F\({dl},{ss}\) = (.*?)(?:   |\n)',proof);assert match
 assert sp.cancel(parse(match[1])-T.forced(dl,ss))==0
 print('Printed F formula PASS',dl,ss,flush=True)
n5=re.search(r'  N_5 = (.*?)\n\n\*\*sl_5',proof,re.S);assert n5
assert sp.expand(parse(n5[1])-T._N5)==0
print('Printed N5 PASS.',flush=True)
source=pickle.load(open('../TOW5_tower_depth10/forced_F910.pkl','rb'))
ours=pickle.load(open('forced_F910.pkl','rb'))
assert set(source)==set(ours)
for key in source:
 if key[0]=='closed': assert sp.sympify(ours[key])==sp.sympify(source[key])
 else: assert sp.cancel(sp.sympify(ours[key])-sp.sympify(source[key]))==0,key
print('Rebuilt forced cache == submitted cache symbolically at every entry.',flush=True)
for dl in range(11):
 assert sp.cancel(T.forced(dl,dl//2)-T2.e_top(dl))==0
 if dl<=6:
  for ss in range(dl//2+1): assert sp.cancel(T.forced(dl,ss)-T3.e(dl,ss))==0
print('Inherited depth<=6 bridge and all tops through 10 PASS.',flush=True)
for dl in (9,10):
 top=sp.factor(pe.ehat(dl,dl//2))
 print('Stable diagonal top',dl,top,flush=True)
 # A nonzero constant or a linear zero below p=4 proves exact degree.
 num,den=sp.fraction(top)
 print('diagonal numerator real roots',sp.solve(num,pe.p),flush=True)
U={dp:sp.cancel(sum((-1)**ss*T.forced(dp,ss) for ss in range(dp//2+1))) for dp in range(11)}
def feed(dl,vs):
 return sum(I*(2*D-1)/(2*D-dp-1)*(T.tau1f(dl-dp)*vs[dp].subs({D:D+1,M:M+1},simultaneous=True)-T.tau2f(dl-dp)*vs[dp].subs({D:D+1,M:M},simultaneous=True)) for dp in range(dl+1))
for dl in (9,10):
 saved=sp.sympify(pickle.load(open(f'a{dl}_closed.pkl','rb')))
 assert sp.cancel(saved-T.arho(dl))==0
 print('Saved a closed form vs Bell PASS',dl,flush=True)
 master=sp.cancel(feed(dl,U)-T.arho(dl));assert master==0
 P=dl//2;last={dp:T.forced(dp,P) for dp in range(dl+1)}
 close=sp.cancel(feed(dl,last)-T.forced(dl,P));assert close==(-1)**P*master==0
 print('Master and closing telescoping identity PASS',dl,flush=True)
for dl in range(2,11):
 aa=T.arho(dl);den=sp.denom(aa)
 numerator=sp.cancel((U[dl]-aa)*den)
 pre=(-1)**(dl//2+1)*2*dl*(dl-1) if dl%2==0 else (-1)**((dl-1)//2)*2*dl*(dl-1)*I*(D-(dl+3)//2)
 Q=sp.cancel(numerator/pre)
 assert sp.Poly(Q,D,M,domain=sp.QQ_I).as_expr()==Q
 assert sp.expand(numerator-pre*Q)==0
 assert sp.degree(Q,D)==3*(dl//2-1)
 print('Defect with shared denominator: exact polynomial Q PASS',dl,'degree',sp.degree(Q,D),flush=True)
# Seed nonuniqueness, not blindness to arbitrary seed errors.
poly=D**5;hs=[]
for ss in range(7):
 hs.append(sp.cancel(poly/T.Pi(5)));poly=sp.expand(poly.subs(D,D+1)-poly)
for ss in range(6):
 H=(2*D-1)/(2*D-11)*hs[ss].subs(D,D+1)
 assert sp.cancel(H-hs[ss]-hs[ss+1])==0
assert hs[6]==0 and sp.cancel(hs[5]-120/T.Pi(5))==0
print('Seed += d^5/Pi5 preserves closing but changes top by 120/Pi5.',flush=True)
h=D**6/T.Pi(5)
for ss in range(6): h=sp.cancel((2*D-1)/(2*D-11)*h.subs(D,D+1)-h)
assert sp.cancel(h-720/T.Pi(5))==0
print('Arbitrary seed blindness is FALSE: seed += d^6/Pi5 changes closing by',sp.factor(h),flush=True)
