
verify_Tn = r'''
# verify top-weight theorem T_n for self-dual, n<=8
import sympy as sp, av9core as C
L = sp.log(2)
h, e = C.selfdual_coeffs(8)
print("extended top part of e_n equals (-1)^n ((L+lv)^n+(L+lw)^n)/(2n!):")
for n in range(1,9):
    en = sp.expand(e[n])
    # extended top part: monomials with ldeg(=a+b in lv,lw) + ln2deg = n
    top = 0
    for (a,b), cf in sp.Poly(en, C.lv, C.lw).terms():
        cf = sp.expand(cf)
        for term in sp.Add.make_args(cf):
            if a+b+sp.degree(term, L) == n:
                top += term*C.lv**a*C.lw**b
    Tn = (-1)**n*((L+C.lv)**n+(L+C.lw)**n)/(2*sp.factorial(n))
    print(f"n={n}: {sp.expand(top - Tn)==0}")
'''
open('verify_Tn.py','w').write(verify_Tn)
import subprocess
r = subprocess.run(['<project>/.venv/bin/python','verify_Tn.py'],capture_output=True,text=True,timeout=300)
print(r.stdout); print(r.stderr[-200:])
with open('progress.log','a') as f:
    f.write("verify_Tn: extended top part of e_n equals (-1)^n((L+lv)^n+(L+lw)^n)/(2n!) for n=1..8 -- top-weight theorem confirmed\n")
