import pickle
import sympy as s
from sympy.polys.fields import field
from sympy.polys.domains import QQ_I
D,M=s.symbols('d m');K,d,m=field((D,M),QQ_I)
raw=pickle.load(open('forced_F910.pkl','rb'))
f={key:K.from_expr(s.sympify(v)) for key,v in raw.items() if isinstance(key[0],int)}
for dl in (9,10):
 a=f[dl,0];u=sum(((-1)**ss*f[dl,ss] for ss in range(dl//2+1)),K.zero)
 assert a.denom==u.denom
 an=a.numer.as_expr();un=u.numer.as_expr()
 degree=s.degree(an,D)
 assert degree=={9:13,10:15}[dl]
 assert s.degree(s.expand(un-an),D)==degree-3
 print('Equal normalized denominators; first three leading d coefficients match:',dl,flush=True)
new=pickle.load(open('N5_proved.pkl','rb'));old=pickle.load(open('../TOW5_tower_depth10/N5_proved.pkl','rb'))
assert s.expand(s.sympify(new['N5'])-s.sympify(old['N5']))==0
print('Freshly regenerated N5 agrees symbolically with submitted N5.',flush=True)
