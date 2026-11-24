import sympy as sp
def logline(msg):
    with open("progress.log","a") as f: f.write(msg+"\n")
s, L, l2 = sp.symbols('s L l2', positive=True)  # L=log pi, l2 = log 2
pi = sp.pi
def Ls(x): return x.subs(sp.log(pi), L).subs(sp.log(2), l2)
ids = {
 'Q(s)  = sum a_n(s)            [f(0) value]': -pi**s/(s+1),
 'P(s)  = sum (-1)^n a_n(s)     [f(pi) value]': pi**s*s/(s+1),
 'R(s)  = sum a_n(s)cos(n pi/2) [f(pi/2)]': (pi/2)**s - pi**s/(s+1),
 'Phi(s)= sum a_n(s)^2          [Parseval f]': 2*pi**(2*s)/(2*s+1) - 2*pi**(2*s)/(s+1)**2,
 'T(s)  = sum b_n(s)sin(n pi/2) [g(pi/2)]': (pi/2)**s,
 'Psi(s)= sum b_n(s)^2          [Parseval g]': 2*pi**(2*s)/(2*s+1),
 'V(s)  = sum b_n(s)/n          [pair w/(pi-x)/2]': pi**(s+1)/((s+1)*(s+2)),
}
out = []
for name, F in ids.items():
    F1 = sp.simplify(sp.diff(F, s)); F2 = sp.simplify(sp.diff(F, s, 2))
    for k in [2,3,4]:
        v  = sp.simplify(Ls(F.subs(s, k)))
        v1 = sp.simplify(Ls(F1.subs(s, k)))
        v2 = sp.simplify(Ls(F2.subs(s, k)))
        out.append((name.split('=')[0].strip(), k, v, v1, v2))
        print(f"{name.split('=')[0].strip():8s} s={k}: F={v}")
        print(f"          F'={v1}")
        print(f"          F''={v2}")
import json
json.dump([[a,k,str(b),str(c),str(d)] for (a,k,b,c,d) in out], open("tables.json","w"))
logline("step_41: exact expansion tables F, F', F'' at s=2,3,4 for Q,P,R,Phi,T,Psi,V")
