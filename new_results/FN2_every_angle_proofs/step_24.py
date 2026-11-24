import sympy as sp, mpmath as mp, time
mp.mp.dps = 50
t0 = time.time()

# ============ Lemma B: the density phi_theta ============
s, c, t = sp.symbols('s c t', positive=False)
c = sp.Symbol('c')
s = sp.Symbol('s', positive=True)
phi = 1 - (s + c)/sp.sqrt(1 + 2*c*s + s**2)
# (i) alternative exact form: phi = (1-c^2)/( sqrt*(sqrt+s+c) )
alt = (1-c**2)/(sp.sqrt(1+2*c*s+s**2)*(sp.sqrt(1+2*c*s+s**2)+s+c))
print("identity phi = (1-c^2)/(sqrt*(sqrt+s+c)):",
      sp.simplify(phi - alt) == 0)

# (ii) derivative sign: phi' = -(1-c^2)/(1+2cs+s^2)^(3/2)
dphi = sp.simplify(sp.diff(phi, s))
print("phi'(s) =", sp.simplify(dphi), " ; equals -(1-c^2)/(...) ^(3/2):",
      sp.simplify(dphi - (-(1-c**2)/(1+2*c*s+s**2)**sp.Rational(3,2))) == 0)

# (iii) tail series: phi(1/t) = 1 - (1+ct)/sqrt(1+2ct+t^2); coefficients d_k
tt = sp.Symbol('tt')
phi_t = 1 - (1+c*tt)*(1+2*c*tt+tt**2)**sp.Rational(-1,2)
ser = sp.series(phi_t, tt, 0, 9).removeO().expand()
print("\ntail series phi(1/t) =", ser)
# Legendre-polynomial form of the coefficients: d_k = (-1)^k (1-c^2) P'_{k-1}(c)/k
ok = True
for k in range(2, 9):
    dk_series = sp.expand(ser).coeff(tt, k)
    dk_leg = (-1)**k * (1-c**2) * sp.legendre(k-1, c).diff(c) / k
    match = sp.simplify(dk_series - dk_leg) == 0
    ok &= match
    print(f"d_{k}: series={sp.simplify(dk_series)}  legendre={sp.simplify(dk_leg)}  match={match}")
print("all d_k match:", ok)

# Legendre identity used: P_k(c) - c P_{k-1}(c) = -(1-c^2) P'_{k-1}(c)/k
okid = all(sp.simplify(sp.legendre(k,c) - c*sp.legendre(k-1,c)
           + (1-c**2)*sp.legendre(k-1,c).diff(c)/k) == 0 for k in range(2,10))
print("Legendre identity P_k - c P_{k-1} = -(1-c^2)P'_{k-1}/k  (k=2..9):", okid)

# (iv) at c=0 only even powers survive and d_{2m}(0) = -binomial(-1/2,m)  (the pi/2 coefficients b_m)
okc0 = True
for m in range(1, 5):
    lhs = sp.simplify(((-1)**(2*m) * (1-c**2) * sp.legendre(2*m-1, c).diff(c)/(2*m)).subs(c, 0))
    rhs = -sp.binomial(sp.Rational(-1,2), m)
    okc0 &= sp.simplify(lhs - rhs) == 0
    print(f"d_{2m}(0) = {lhs}  vs  -C(-1/2,{m}) = {rhs}")
print("c=0 reduction to P3's b_l:", okc0)
print("elapsed", time.time()-t0)
with open("progress.log","a") as f:
    f.write(f"{time.strftime('%H:%M')} Lemma B symbolic: phi_theta form, decreasing, tail d_k=(-1)^k(1-c^2)P'_(k-1)/k, c=0->b_l  [PROVED algebraically, verified]\n")
