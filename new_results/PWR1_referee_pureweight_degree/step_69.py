import mpmath as mp, pickle, sympy as sp
mp.mp.dps = 105
exec(open('step_67.py').read().split("# K3's stored Qpolys")[0])  # reuse defs (Dtilde, Theta)
Q = pickle.load(open('../K3_fourier_all_orders/Qpolys.pkl','rb'))
u_sym = sp.Symbol('u')
zsym = {f'z{k}': sp.Symbol(f'z{k}') for k in range(2, 13)}
znum = {sp.Symbol(f'z{k}'): mp.zeta(k) for k in range(2, 13)}
print("Theta_n(u) [my D~-build] vs K3's PSLQ-identified Q_n(u), at u = 0,1,2,5,7:")
for n in range(2, 9):
    qn = sp.sympify(Q[n], locals=zsym)
    diffs = []
    for uv in [0, 1, 2, 5, 7]:
        a = Theta(n, uv)
        b = mp.mpf(str(qn.subs(u_sym, uv).subs(znum).evalf(50))) if not isinstance(qn, sp.Symbol) else 0
        b = mp.mpf(str(qn.subs({u_sym: uv, **{sp.Symbol(f'z{k}'): sp.zeta(k) for k in range(2,13)}}).evalf(50)))
        diffs.append(abs(a-b))
    print(f"  n={n}: max |diff| over u in {{0,1,2,5,7}} = {mp.nstr(max(diffs), 3)}")
# leading-coefficient identity check (symbolic): 2(-4)^{m-1}(1-2^{1-2m}) =?= (-1)^{m-1}(2^{2m-1}-1)
m = sp.symbols('m', integer=True, positive=True)
for mv in range(1, 9):
    lhs = 2*(-4)**(mv-1)*(1-sp.Rational(2)**(1-2*mv))
    rhs = (-1)**(mv-1)*(2**(2*mv-1)-1)
    assert sp.simplify(lhs-rhs) == 0, mv
print("leading-coefficient identity 2(-4)^{m-1}(1-2^{1-2m}) = (-1)^{m-1}(2^{2m-1}-1): TRUE for m=1..8 (exact)")
# and Theta top coefficient vs paper's claim for n=2m <= 12:
print("\npaper claim: leading coeff of Theta_{2m} = (-1)^{m-1}(2^{2m-1}-1) zeta(2m):")
for mv in range(1, 7):
    n = 2*mv
    top = Dtilde(n, mv-1)*(-4)**(mv-1)/mp.factorial(2*mv-2)
    claim = (-1)**(mv-1)*(2**(2*mv-1)-1)*mp.zeta(2*mv)
    print(f"  m={mv}: {mp.nstr(top, 30)}  diff={mp.nstr(abs(top-claim),3)}")
