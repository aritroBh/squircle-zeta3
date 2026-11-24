exec(open('common.py').read())
forms = load_forms()
C4, Q7 = get_C4_Q7(forms)
P8 = dict_to_sym(forms[8][3])
print("C4 quotient exact:", sp.expand(a2*a3*a4*SA*C4/10 - P8)==0)
Pc4 = sp.Poly(C4, a1,a2,a3,a4,a5)
print("C4 number of terms:", len(Pc4.terms()))
# z2z(3,5): check ratio
P235 = dict_to_sym(forms[10][5])
# candidate factors
cand1 = -a2*a3*a4*SA*(a2+a4)*(a1+a3+a5)*C4/10
print("identity as stated:", sp.expand(P235 - cand1)==0)
# Try without minus:
cand2 = a2*a3*a4*SA*(a2+a4)*(a1+a3+a5)*C4/10
print("identity without minus:", sp.expand(P235 - cand2)==0)
# Try factoring P235 directly
f235 = sp.factor(P235)
print("\nfactor(P235):", f235)
