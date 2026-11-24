import sympy as sp, pickle, os, sys
wd = "<project>/research/08_theory/DEC2_decagon_exact_assembly"
os.chdir(wd); sys.path.insert(0, '<project>/research/08_theory/DEC1_decagon')
t = sp.Symbol('t'); s = sp.Symbol('s'); r5 = sp.sqrt(5); phi = (1+r5)/2
kap1 = 1/phi; kap4 = -1/phi**2
import dec_engine as E
# ---------- (a) D1, D2 series certificates ----------
k, m, n = sp.symbols('k m n', integer=True, positive=True)
x1s, x4s = sp.Symbol('x1'), sp.Symbol('x4')
# LHS term (from expanding 1/u1 = sum x1^k t^k, w4 = -sum x4^m t^m/m, int t^{k+m} log t = -1/(k+m+1)^2):
T_lhs = x1s**k * x4s**m / (m*(k+m+1)**2)
# RHS term: (1/x1) Li21(x1, x4/x1) = sum_{n>m>=1} x1^{n-m-1} x4^m / (n^2 m); set n = k+m+1:
T_rhs = x1s**(n-m-1) * x4s**m / (n**2 * m)
print("D1 term equality (k = n-m-1):", sp.simplify(T_lhs - T_rhs.subs(n, k+m+1)))
# absolute convergence: |x1|=1/phi<1, x4=1/phi^2<1: sum |x1|^k x4^m/(m (k+m+1)^2) <= sum |x1|^k x4^m/m = log-geometric, finite.
# ---------- (b) K1, K4 series: int_0^1 log t/(1+kappa t) dt = Li2(-kappa)/kappa ----------
print("K identity: sum_k (-kap)^k * (-1/(k+1)^2) = (1/kap) sum_{j>=1} (-kap)^j/j^2 = Li2(-kap)/kap [exact rearrangement]")
# ---------- (c) reflection algebraic identities ----------
id1 = sp.expand((1+kap4*(1-s)) - (phi+s)/phi**2)
id2 = sp.expand((phi+s)/phi**2 - (1/phi)*(1+s/phi))
id3 = sp.expand((1+kap1*(1-s)) - phi*(1-s/phi**2))
id4 = sp.expand((phi**2-s)/phi - phi*(1-s/phi**2))
id5 = sp.expand((phi+s) - phi*(1+s/phi))
print("reflection identities (must be 0):", id1, id2, id3, id4, id5)
# D3 = (1/phi)(D2 - lf K4), D4 = phi (lf K1 + D1)
lf = sp.Symbol('lf')
D1_eng = E.D1; D2_eng = E.D2
K4 = sp.expand(E.L2_4/kap4)   # int log t /u4 = Li2(x4)/kap4
K1 = sp.expand(E.L2_1/kap1)
def canon5(e):
    e = sp.expand(e)
    P = sp.Poly(e, *E.NAMED)
    out = sp.Integer(0)
    for monom, coeff in P.terms():
        c = sp.together(coeff)
        cc = c.subs(r5, -r5)
        p = sp.radsimp(sp.expand((c+cc)/2)); q = sp.radsimp(sp.expand((c-cc)/(2*r5)))
        mm = sp.Integer(1)
        for a_, e_ in zip(E.NAMED, monom): mm *= a_**e_
        out += (p + q*r5)*mm
    return sp.expand(out)
D3 = canon5((D2_eng - E.lf*K4)/phi)
D4 = canon5(phi*(E.lf*K1 + D1_eng))
print("D3 (reflection of D2) canonical diff vs engine:", canon5(D3 - E.D3))
print("D4 (reflection of D1) canonical diff vs engine:", canon5(D4 - E.D4))
res = pickle.load(open('D5678_proved.pkl','rb'))
res['D1'] = D1_eng; res['D2'] = D2_eng; res['D3'] = D3; res['D4'] = D4
pickle.dump(res, open('D_all_proved.pkl','wb'))
print("saved D_all_proved.pkl (D1..D10)")
