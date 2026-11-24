
import mpmath as mp
mp.mp.dps = 60
sq2 = mp.sqrt(2)
a = 1+sq2; b = sq2-1
z3 = mp.zeta(3); pi = mp.pi; lam = mp.log(2)

def t_of_v(v): return (v*v+2*v-1)/(4*v)
def omt_of_v(v): return (1+2*v-v*v)/(4*v)

def R20(k,v): return (1-sq2)/(4*v) + 1/(4*(v-a)) + (2*k-7)/(4*(v+a)) - (k-2)/(2*(v+b))
def R11(k,v): return (2*k-5+sq2)/(4*v) - (2*k-5)/(2*(v+a))
def S20(v): return 1/(2*(v+a)) - 1/(2*(v+b))
def S11(v): return 1/(2*v) - 1/(v+a)

def W3(k):
    f = lambda v: R20(k,v)*mp.log(t_of_v(v))**2 + R11(k,v)*mp.log(t_of_v(v))*mp.log(omt_of_v(v))
    return mp.quad(f, [b, a])
def slope_int():
    f = lambda v: S20(v)*mp.log(t_of_v(v))**2 + S11(v)*mp.log(t_of_v(v))*mp.log(omt_of_v(v))
    return mp.quad(f, [b, a])

slope = slope_int()
target_slope = -(z3/8 + pi**2*lam/24 + lam**3/12)
print("slope integral      :", mp.nstr(slope, 40))
print("target -(z3/8+...)  :", mp.nstr(target_slope, 40))
print("residual            :", mp.nstr(slope-target_slope, 6))

W = {k: W3(k) for k in [3,15,16,18,20]}
for k in sorted(W): print(f"W3({k}) =", mp.nstr(W[k],30))
# affine law: W3(k) - W3(3) = (k-3)*slope
W3p = W[3]
for k in [15,16,18,20]:
    lhs = W[k]-W3p; rhs=(k-3)*slope
    print(f"k={k}: (W3(k)-W3(3))-(k-3)*slope =", mp.nstr(lhs-rhs,6))
