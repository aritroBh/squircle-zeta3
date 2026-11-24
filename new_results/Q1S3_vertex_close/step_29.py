import importlib, q3_engine, q1s_exact
importlib.reload(q1s_exact); importlib.reload(q3_engine)
from q3_engine import *
from q1s_exact import *
import mpmath as mp
mp.mp.dps = 60
def fembc(x): return femb(x,1)
b = mp.sqrt(2)-1; a=1+mp.sqrt(2)

# Control 1 target: int P/(1-by)^2, P = log y * log(1-by)
val_direct = mp.quad(lambda y: mp.log(y)*mp.log(1-b*y)/(1-b*y)**2, [0,1])
print("direct:", mp.nstr(val_direct,20))
# term A: int a^2 log(1-y/a) / (y(y-a)) dy
tA = mp.quad(lambda y: a*a*mp.log(1-y/a)/(y*(y-a)), [0,1])
print("termA direct:", mp.nstr(tA,20))
# term B: int a^2 log y/(y-a)^2 dy
tB = mp.quad(lambda y: a*a*mp.log(y)/(y-a)**2, [0,1])
print("termB direct:", mp.nstr(tB,20), " a log(1-1/a) =", mp.nstr(a*mp.log(1-1/a),20))
print("sum:", mp.nstr(tA+tB,20))
# so the split: int P/(1-by)^2 = -(1/b) int P'/(1-by) dy; P' = (1/y)log(1-by) + log y * (-b/(1-by))
# -(1/b)*(-a) int [..]/(y-a) = a^2 int [(1/y)log(1-y/a) + log y/(y-a)]/(y-a) dy = tA + tB  -- wait!
# term B integrand should be a^2 * log y / ((y-a)(y-a)) -- that is (1/(y-a))*(1/(y-a)): sF=a,sK=a square OK
# term A: a^2 log(1-y/a)/(y(y-a)): sF=0,sK=a
print("tA+tB - direct:", mp.nstr(tA+tB-val_direct,10))
