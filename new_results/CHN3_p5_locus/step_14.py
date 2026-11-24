exec(open('common.py').read())
forms = load_forms()
C4, Q7 = get_C4_Q7(forms)
print("C4: terms", len(sp.Poly(C4).terms()), " C4(1,1,1,1,1) =", C4.subs({v:1 for v in A}))
print("Q7': terms", len(sp.Poly(Q7).terms()))
# verify z2z(3,5) identity
P235 = dict_to_sym(forms[10][5])
print("z2z(3,5) == -a2a3a4(SA)(a2+a4)(a1+a3+a5)C4/10 :",
      sp.expand(P235 + a2*a3*a4*SA*(a2+a4)*(a1+a3+a5)*C4/10) == 0)
# Save C4, Q7 as strings for reuse
open('C4.txt','w').write(str(C4))
open('Q7p.txt','w').write(str(Q7))
# Restrictions: C4 on a1=0 should be P4's C3 (in a2..a5) times ... let me check what CHN1's C3 is
# From CHN1 result: P8^(4) = a1a2a3a4(sum a) C3/10 in vars a1..a4. C4|a1=0 = ? in vars a2..a5.
print("\nC4|a1=0 factored:", sp.factor(C4.subs(a1,0)))
print("\nC4|a5=0 factored:", sp.factor(C4.subs(a5,0)))
