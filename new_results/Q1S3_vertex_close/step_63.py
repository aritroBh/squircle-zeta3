import importlib, q3_engine, q1s_exact, q3_eval
for m in (q1s_exact, q3_engine, q3_eval): importlib.reload(m)
from q3_engine import *
from q3_eval import GEval
from q1s_exact import femb
import mpmath as mp, time
mp.mp.dps=80
def fembc(x): return femb(x,1)
r2=mp.sqrt(2); a=1+r2; b=r2-1
ge=GEval(80, M=48, gln=4, subn=4)
t0=time.time()
tot = ge.rest_total(fembc(smb), an if (an:=fembc(sa)) else None)
print("rest(-b,a) total:", mp.nstr(tot,18), " time:", round(time.time()-t0,1))
# exact via series-ish: G(-b,a;1) = int_0^1 log(1-t/a)/(t+b) dt
ex = mp.quad(lambda t: mp.log(1-t/a)/(t+b), [0,1])
print("direct:", mp.nstr(ex,18), " resid:", mp.nstr(abs(tot-ex),4))
arr=ge.restF[(fembc(smb),fembc(sa))]
# F at x=0.31:
j=min(range(len(ge.pts)), key=lambda j: abs(ge.pts[j][0]-mp.mpf('0.31')))
x=ge.pts[j][0]
ex2=mp.quad(lambda t: mp.log(1-t/a)/(t+b), [0,x])
print("F(x=%s):"%mp.nstr(x,4), mp.nstr(arr[j],18), mp.nstr(ex2,18), " resid:", mp.nstr(abs(arr[j]-ex2),4))
