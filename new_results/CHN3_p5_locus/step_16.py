exec(open('common.py').read())
forms = load_forms()
P8 = dict_to_sym(forms[8][3])
P10 = dict_to_sym(forms[10][6])
P235 = dict_to_sym(forms[10][5])

def exact_div(num, den):
    q, r = sp.div(num, den, a1,a2,a3,a4,a5)
    assert sp.expand(r)==0, ("remainder nonzero!", r)
    return sp.expand(q)

C4 = exact_div(10*P8, a2*a3*a4*SA)
Q7 = exact_div(224*P10, a2*a3*a4*SA)
print("C4 terms:", len(sp.Poly(C4,a1,a2,a3,a4,a5).terms()), " C4(1..1) =", C4.subs({v:1 for v in A}))
print("Q7' terms:", len(sp.Poly(Q7,a1,a2,a3,a4,a5).terms()))

# z2z(3,5) identity
rem1 = sp.expand(P235 + a2*a3*a4*SA*(a2+a4)*(a1+a3+a5)*C4/10)
print("identity as stated (with -):", rem1==0)
rem2 = sp.expand(P235 - a2*a3*a4*SA*(a2+a4)*(a1+a3+a5)*C4/10)
print("identity with +:", rem2==0)

open('C4.txt','w').write(str(C4))
open('Q7p.txt','w').write(str(Q7))
