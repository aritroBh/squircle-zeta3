import sympy as sp
z3, z5, z7, pi_s = sp.symbols('z3 z5 z7 pi')
def sym_zeta(k):
    return {3:z3,5:z5,7:z7}.get(k, sp.zeta(k))
def sym_eta(k): return (1-2**(1-k))*sym_zeta(k)
def terms_a(s):
    if s==1: return [(2/pi_s,'a',2),(-2/pi_s,'c',2)]
    if s==2: return [(sp.Integer(4),'a',2)]
    if s==3: return [(6*pi_s,'a',2),(12/pi_s,'c',4),(-12/pi_s,'a',4)]
    if s==4: return [(8*pi_s**2,'a',2),(sp.Integer(-48),'a',4)]
def terms_b(s):
    if s==1: return [(sp.Integer(-2),'a',1)]
    if s==2: return [(-2*pi_s,'a',1),(-4/pi_s,'c',3),(4/pi_s,'a',3)]
    if s==3: return [(-2*pi_s**2,'a',1),(sp.Integer(12),'a',3)]
    if s==4: return [(-2*pi_s**3,'a',1),(24*pi_s,'a',3),(48/pi_s,'c',5),(-48/pi_s,'a',5)]
def Ssum(w,k): return sym_zeta(k) if w=='c' else -sym_eta(k)
def B_sym(s,q):
    tot = sp.Integer(0)
    for (c1,w1,k1) in terms_a(s):
        for (c2,w2,k2) in terms_b(q):
            w = 'c' if w1==w2 else 'a'
            tot += sp.sympify(c1)*sp.sympify(c2)*Ssum(w,k1+k2)
    return sp.simplify(tot)
out = {}
for (s,q) in [(1,1),(1,2),(2,1),(2,2),(1,3),(2,3),(3,3),(1,4),(2,4),(3,4),(4,4)]:
    expr = sp.simplify(B_sym(s,q))
    out[(s,q)] = expr
    print(f"B({s},{q}) = {expr}")
import json
json.dump({f"{s},{q}": str(e) for (s,q),e in out.items()}, open("B_forms.json","w"))
