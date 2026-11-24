import sympy as sp
u = sp.symbols('u'); k = sp.symbols('kappa', positive=True)
S = sp.sqrt(u**2 + 4*k**2); D = sp.sqrt(1+4*k**2)
mu_w1   = sp.Rational(1,2)*(1/S + 1/(u-1) + D/(S*(u-1)))
mu_wm1  = sp.Rational(1,2)*(1/S + 1/(u+1) - D/(S*(u+1)))
mu_wr   = sp.Rational(1,2)*(1/S + 1/(u-1) - D/(S*(u-1)))
mu_w0   = 1/S
# corner measures (same combinations as LEAD2)
m_odd  = sp.simplify(mu_wm1 - mu_w1)          # j < (d-1)/2 odd d
m_e    = sp.simplify(mu_wm1 - mu_wr)          # even d
m_diag = sp.simplify(mu_wm1 - mu_w0/2)        # diagonal
m_oddlast = sp.simplify(mu_wm1 - mu_w1/2 + mu_wr/2 - mu_w0/2)  # oddlast: 1/(v+a)-1/(2(v-a))+1/(2(v+b))-1/(2v) in w: mu_{-1} - mu_1/2 + mu_{-r}/2 - mu_0/2
print("m_odd  =", sp.simplify(m_odd.rewrite(sp.sqrt)))      # expect (1 + D u/S)/(1-u^2)
print("target :", (1 + D*u/S)/(1-u**2))
print("m_e    =", sp.simplify(m_e))                         # expect (1 - D/S)/(1-u^2)
print("target :", (1 - D/S)/(1-u**2))
print("m_diag =", sp.simplify(m_diag))                      # expect (1-u)/2 m_e
print("target :", sp.simplify((1-u)/2*(1 - D/S)/(1-u**2)))
print("m_oddlast =", sp.simplify(m_oddlast))
print("target     :", sp.simplify(m_odd - 1/(2*(1-u))))
