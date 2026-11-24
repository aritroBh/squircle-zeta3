
import pickle, time
import sympy as sp
t = sp.Symbol('t')
el0, el1 = sp.symbols('el0 el1')       # log t, log(1-t)
W1, W4 = sp.symbols('W1 W4')           # log(1+t/phi), log(1-t/phi^2)
L2, LF = sp.symbols('L2 LF')           # log 2, log phi
sh, hh = sp.symbols('sh hh')
l0, l1 = sp.symbols('l0 l1')
taus = sp.symbols('tau1:6'); gams = sp.symbols('gam1:6')
al, be = sp.symbols('alpha beta')
loc = {'t':t,'alpha':al,'beta':be,'l0':l0,'l1':l1,'sh':sh,'hh':hh}
for i,x in enumerate(taus): loc[f'tau{i+1}']=x
for i,x in enumerate(gams): loc[f'gam{i+1}']=x
E5oL = sp.sympify(pickle.load(open("E5_dec_stage1.pkl","rb")), locals=loc)

s5 = sp.sqrt(5)
phi = (1+s5)/2
B1 = (3-s5)/4        # 1/(2 phi^2)
B4 = (s5-1)/4        # 1/(2 phi)
ap = (s5-2)/2        # 1/(2 phi^3)
lam1 = W1 - L2 - 2*LF
lam4 = W4 - L2 - LF
lam0 = -L2 - LF
a1t = B1*(1+t/phi); a4t = B4*(1-t/phi**2); a0c = B4
subsdict = {l0: el0 - L2 - 2*LF, l1: el1 - L2 - 2*LF,
            al: B1, be: B1, sh: -sp.Rational(1,2), hh: (5+2*s5)/4}
for k in range(1,6):
    subsdict[taus[k-1]] = a1t*lam1**k + a4t*lam4**k + a0c*lam0**k
    subsdict[gams[k-1]] = ap*lam1**k - ap*lam4**k
t0=time.time()
E5d = E5oL.subs(subsdict)
print("subs %.1fs"%(time.time()-t0))
pickle.dump(str(E5d), open("E5_dec_subs.pkl","wb"))
print("saved E5_dec_subs.pkl, len", len(str(E5d)))
