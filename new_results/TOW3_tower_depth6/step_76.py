import sympy as sp
s2 = sp.sqrt(2); ii = sp.I
aq = 1+s2; bq = s2-1
def w_parts(mm):
    w = sp.expand((-aq-ii)**(-mm) - (bq-ii)**(-mm))
    w = sp.radsimp(w)
    wre = sp.expand(sp.re(w)); wim = sp.expand(sp.im(w))
    wr0 = sp.Rational(wre.coeff(s2,0)); ws0 = sp.Rational(wre.coeff(s2,1))
    wr1 = sp.Rational(wim.coeff(s2,0)); ws1 = sp.Rational(wim.coeff(s2,1))
    return (wr0,wr1),(ws0,ws1)
for mm in range(1,5):
    print(mm, w_parts(mm))
