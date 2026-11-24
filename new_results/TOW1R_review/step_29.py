
import mpmath as mp
mp.mp.dps = 80
sq2=mp.sqrt(2); a=1+sq2; b=sq2-1
z3=mp.zeta(3); pi=mp.pi; lam=mp.log(2)
def t(v): return (v*v+2*v-1)/(4*v)
def om(v): return (1+2*v-v*v)/(4*v)
def R20(k,v): return (1-sq2)/(4*v)+1/(4*(v-a))+(2*k-7)/(4*(v+a))-(k-2)/(2*(v+b))
def R11(k,v): return (2*k-5+sq2)/(4*v)-(2*k-5)/(2*(v+a))
def S20(v): return 1/(2*(v+a))-1/(2*(v+b))
def S11(v): return 1/(2*v)-1/(v+a)
def W3(k): return mp.quad(lambda v: R20(k,v)*mp.log(t(v))**2+R11(k,v)*mp.log(t(v))*mp.log(om(v)),[b,a])
slope = mp.quad(lambda v: S20(v)*mp.log(t(v))**2+S11(v)*mp.log(t(v))*mp.log(om(v)),[b,a])
tgt = -(z3/8+pi**2*lam/24+lam**3/12)
print("slope residual (dps80):", mp.nstr(slope-tgt,5))
W={k:W3(k) for k in [3,17,19]}
for k in [17,19]:
    print(f"k={k}: (W3(k)-W3(3))-(k-3)*slope =", mp.nstr(W[k]-W[3]-(k-3)*slope,5))
