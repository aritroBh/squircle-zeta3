import sympy as sp
s1,s2,s3,s4,s5=sp.symbols('s1:6')
S=(s1,s2,s3,s4,s5)
s13=s4-s1-s2; s24=s5-s2-s3
P2_11=-s3*s4+s1*(s3-s5)
m11=s3*(-s1*(s1+2*s2+s3)+s3*s4+s4**2)+s1*s5*(s1+s5)
m12=-s13*s24*(s1+s2+s3+s4+s5)
m21=s1*s3*(s1+s2+s3-2*(s4+s5))
def shA(f): return f.subs({s1:s2,s2:s3,s3:s4,s4:s5,s5:s1},simultaneous=True)
def shB(f): return f.subs({s1:s5,s2:s1,s3:s2,s4:s3,s5:s4},simultaneous=True)
for name,sh in (('A: (s2,s3,s4,s5,s1)',shA),('B: (s5,s1,s2,s3,s4)',shB)):
    for lab,f,target in (('P2',P2_11,s13*s24),('M3',m11,m12)):
        cand=sp.cancel(s13/(s4-s1)*(f-sh(f)))
        print(name,lab,'matches target:',sp.simplify(cand-target)==0)
# relabeling 2<->3 consistency: m21 = m12 relabeled
rel={s1:s13,s3:s24}  # s2,s4,s5 fixed
print('m21 == m12|2<->3 :',sp.simplify(m12.subs({s1:s13,s3:s24},simultaneous=True)-m21)==0)
