"""indep_L.py (referee): derive L_r = [x^r] log A(i+x), r<=8, DIRECTLY by series expansion (no Faa di Bruno, no tow2lib.Lr),
from A(v) = C v^(2d-1) S_m(t(v)) (v+i)^(2-2d) / (v^4-6v^2+1),  t(v)=v/4+1/2-1/(4v)  (CX9's form).
Only the PROVED S-derivative ratios R_1..R_4 (via tow2lib.sl cumulants) are imported; compare with L_simplified cache."""
import sys, pickle, time
sys.path.insert(0, 'lib')
import sympy as sp
import tow2lib as T2
d, m = T2.d, T2.m
x = sp.Symbol('x')
I = sp.I
ORD = 8
def trunc(expr, n=ORD):
    return sp.series(sp.expand(expr), x, 0, n+1).removeO()
t0 = time.time()
# algebraic pieces
P = sp.expand((I+x)**4 - 6*(I+x)**2 + 1)
assert sp.simplify(P.subs(x,0)) == 8
y = sp.expand(P/8 - 1)
logP = sum((-1)**(n-1) * y**n / n for n in range(1, ORD+1))
logP = trunc(logP)
log1 = sum((-1)**(n-1) * (x/I)**n / n for n in range(1, ORD+1))        # log(1 + x/i)
log2 = sum((-1)**(n-1) * (x/(2*I))**n / n for n in range(1, ORD+1))    # log(1 + x/(2i))
# h(x) = t(i+x) - t(i)
hx = trunc(sp.series(x/4 - sp.Rational(1,4)*(1/(I+x) - 1/I), x, 0, ORD+1).removeO())
assert sp.expand(hx).coeff(x,1) == 0 and sp.expand(hx).coeff(x,0) == 0   # t'(i)=0
logS = 0
for j in range(1, 5):
    logS += T2.sl[j] * trunc(hx**j) / sp.factorial(j)
logA = (2*d-1)*log1 + (2-2*d)*log2 - logP + logS
logA = sp.expand(logA)
Lraw = pickle.load(open('lib/L_simplified.pkl','rb'))
print('series built', round(time.time()-t0,1), 's', flush=True)
allok = True
for r in range(1, ORD+1):
    ind = logA.coeff(x, r)
    cache = sp.sympify(Lraw[r], locals={'d':d,'m':m,'I':sp.I})
    diff = sp.cancel(sp.together(ind - cache))
    ok = (diff == 0); allok &= ok
    print(f'r={r}: independent series L_r == L_simplified[r]: {ok}', flush=True)
print('ALL r=1..8 equal:', allok, round(time.time()-t0,1), 's')
