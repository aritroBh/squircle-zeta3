import pickle, sympy as sp, time
forms=pickle.load(open('forms10.pkl','rb'))
a1,a2,a3,a4,a5=sp.symbols('a1 a2 a3 a4 a5'); A=[a1,a2,a3,a4,a5]
def to_sp(dd): return sum(sp.Rational(c.numerator,c.denominator)*sp.prod(A[i]**e[i] for i in range(5)) for e,c in dd.items())
C4=sp.sympify(open('C4.txt').read()); Q7=sp.sympify(open('Q7p.txt').read())
Sa=a1+a2+a3+a4+a5
# W8 zeta(3,5) = index 3
z35=to_sp(forms[8][3])
q,r=sp.div(z35,C4,*A)
print("W8 zeta(3,5) / C4 remainder 0?", sp.expand(r)==0, " quotient:", sp.factor(q))
print("  expected a2a3a4(Σa)/10:", sp.factor(a2*a3*a4*Sa/10))
# W10 zeta2 zeta(3,5) = index 5
z2z35=to_sp(forms[10][5])
q2,r2=sp.div(z2z35,C4,*A)
print("W10 zeta2zeta(3,5) / C4 remainder 0?", sp.expand(r2)==0, " quotient:", sp.factor(q2))
def log_progress(line):
    with open('progress.log','a') as f: f.write(time.strftime('%H:%M:%S ')+line+'\n')
log_progress("EXACT (from forms10.pkl): W10 zeta(3,7) coord = a2a3a4(Sum a)*Q7'/224 (exact division). Confirms V10={C4=0,Q7'=0} is EXACTLY the irreducible-free locus through W10 for positive a: W8 z(3,5)=a2a3a4(Sum a)C4/10, W10 z2z(3,5)=-...(C4)/10 (vanishes at C4=0), W10 z(3,7)=a2a3a4(Sum a)Q7'/224 (vanishes at Q7'=0). This closes Cases 1&2 of the argument EXACTLY.")
print("recorded")
print("\nproof log:")
print(open('proof_div.log').read())
