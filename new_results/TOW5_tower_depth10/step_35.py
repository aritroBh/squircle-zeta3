import sympy as sp, time, pickle
t0=time.time()
work = '<project>/research/08_theory/TOW5_tower_depth10'
d, m = sp.symbols('d m')
i = sp.I
raw = pickle.load(open(work+'/R5_raw.pkl','rb'))
num5 = sp.sympify(raw['num5'], locals={'d':d,'m':m,'I':i})
den5 = sp.sympify(raw['den5'], locals={'d':d,'m':m,'I':i})
R5 = num5/den5

# R_1..R_4 from tow2lib
import sys; sys.path.insert(0, work)
import tow2lib as T2
R = {l: T2.R[l] for l in range(1,5)}

# cumulant/moment relation: verify formula for sl_5 via series reversion
h = sp.Symbol('h')
Rser = {l: sp.Symbol(f'R{l}') for l in range(1,6)}
f = 1 + sum(Rser[l]*h**l/sp.factorial(l) for l in range(1,6))
logf = sp.series(sp.log(f), h, 0, 6).removeO().expand()
sl5_sym = sp.expand(logf.coeff(h,5)*sp.factorial(5))
sl4_sym = sp.expand(logf.coeff(h,4)*sp.factorial(4))
print("sl_5 formula:", sl5_sym)
print("sl_4 formula matches tow2lib:", sp.expand(sl4_sym - (Rser[4]-4*Rser[1]*Rser[3]-3*Rser[2]**2+12*Rser[1]**2*Rser[2]-6*Rser[1]**4))==0)

sl5 = sp.together(sl5_sym.subs({Rser[l]: (R5 if l==5 else R[l]) for l in range(1,6)}))
sn, sd = sp.fraction(sl5)
print("sl_5 denominator factored:", sp.factor(sd))
print("sl_5 numerator degree in d:", sp.degree(sp.expand(sn), d), f"({round(time.time()-t0,1)}s)")
