import sympy as sp
i = sp.I
sq2 = sp.sqrt(2)
a = 1+sq2; b = sq2-1
def tau1(n):
    s=n+1
    return sp.simplify((-1)**s * ( (i-b)**(-(s+1)) + (i+a)**(-(s+1)) - i**(-(s+1)) ))
def tau2(n):
    s=n+1
    return sp.simplify((-1)**s * ( (i-a)**(-(s+1)) + (i+b)**(-(s+1)) - i**(-(s+1)) ))
# proof_common dict values
tdict = { (1,0): (-1-i)/2, (2,0): (-1+i)/2, (1,1): sp.Rational(3,4)*(1-i),
      (2,1): -sp.Rational(3,4)-sp.Rational(3,4)*i, (1,2): 1+sp.Rational(3,4)*i, (2,2): 1-sp.Rational(3,4)*i }
print("tau1_0 =", tau1(0), " match:", sp.simplify(tau1(0)-tdict[(1,0)])==0)
print("tau1_1 =", tau1(1), " match:", sp.simplify(tau1(1)-tdict[(1,1)])==0)
print("tau1_2 =", tau1(2), " match:", sp.simplify(tau1(2)-tdict[(1,2)])==0)
print("tau2_0 =", tau2(0), " match:", sp.simplify(tau2(0)-tdict[(2,0)])==0)
print("tau2_1 =", tau2(1), " match:", sp.simplify(tau2(1)-tdict[(2,1)])==0)
print("tau2_2 =", tau2(2), " match:", sp.simplify(tau2(2)-tdict[(2,2)])==0)
# check tc1c: conj(tau1_n) = (-1)^n tau2_n
for n in range(4):
    lhs = sp.conjugate(tau1(n)); rhs = (-1)**n * tau2(n)
    print(f"conj(tau1_{n}) == (-1)^{n} tau2_{n}:", sp.simplify(lhs-rhs)==0, "  (tau1_%d=%s, tau2_%d=%s)"%(n,tau1(n),n,tau2(n)))
# check tau2 = conj(tau1)?
for n in range(3):
    print(f"tau2_{n} == conj(tau1_{n}):", sp.simplify(tau2(n)-sp.conjugate(tau1(n)))==0)
# verify DL1(v) = sum tau1_n (v-i)^{n+1}: DL1 = 1/(v-b)+1/(v+a)-1/v, DL1(i)=0?
v=sp.symbols('v')
DL1 = 1/(v-b)+1/(v+a)-1/v
DL2 = 1/(v-a)+1/(v+b)-1/v
print("DL1(i)=", sp.simplify(DL1.subs(v,i)), " DL2(i)=", sp.simplify(DL2.subs(v,i)))
# tau1_0 should be DL1'(i)
print("DL1'(i)=", sp.simplify(sp.diff(DL1,v).subs(v,i)), "== tau1_0:", sp.simplify(sp.diff(DL1,v).subs(v,i)-tau1(0))==0)
print("DL2'(i)=", sp.simplify(sp.diff(DL2,v).subs(v,i)), "== tau2_0:", sp.simplify(sp.diff(DL2,v).subs(v,i)-tau2(0))==0)
# tau1_1 = DL1''(i)/2
print("DL1''(i)/2==tau1_1:", sp.simplify(sp.diff(DL1,v,2).subs(v,i)/2 - tau1(1))==0)
print("DL1'''(i)/6==tau1_2:", sp.simplify(sp.diff(DL1,v,3).subs(v,i)/6 - tau1(2))==0)
with open("<project>/research/08_theory/TOW1R_review/progress.log","a") as fh:
    fh.write("step_06: verified tau1_n,tau2_n formula + proof_common t-dict values + tc1c conj relation\n")
