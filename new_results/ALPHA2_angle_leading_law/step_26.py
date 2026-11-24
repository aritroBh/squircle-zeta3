import sympy as sp
u = sp.symbols('u'); k = sp.symbols('kappa', positive=True)
S = sp.sqrt(u**2 + 4*k**2); D = sp.sqrt(1+4*k**2)
m_odd = (1 + D*u/S)/(1-u**2); m_e = (1 - D/S)/(1-u**2)
# key LEAD2 identities for the corner measures
print("m_odd - (1/(1-u) - u m_e) =", sp.simplify(m_odd - (1/(1-u) - u*m_e)))
m_diag = sp.simplify((1-u)/2*m_e)
print("m_diag - m_e/2 + u m_e/2 =", sp.simplify(m_diag - m_e/2 + u*m_e/2))
# m_oddlast - (m_odd - 1/(2(1-u)))
m_oddlast = sp.simplify(m_odd - 1/(2*(1-u)))
print("ok")
