src = open('ftfield.py').read()
src = src.replace("    C = A1*C2+C1*A2 - (B1*C2+C1*B2) + 3*C1*C2\n    return (A,B,C,d1*d2)",
                  "    C = A1*C2+C1*A2 + B1*B2 - (B1*C2+C1*B2) + 3*C1*C2\n    return (A,B,C,d1*d2)")
# lazy high-precision embeddings
src = src.replace("""import mpmath as mp
EMB = [2*mp.cos(2*mp.pi/7), 2*mp.cos(4*mp.pi/7), 2*mp.cos(6*mp.pi/7)]""",
"""import mpmath as mp
_EMB = None
def get_EMB(prec=60):
    global _EMB
    old = mp.mp.dps; mp.mp.dps = prec
    _EMB = [2*mp.cos(2*mp.pi/7), 2*mp.cos(4*mp.pi/7), 2*mp.cos(6*mp.pi/7)]
    mp.mp.dps = old
    return _EMB
EMB = get_EMB()""")
open('ftfield.py','w').write(src)
test = r'''
import ftfield as F
a = F.AL
z3 = F.fsub(F.fpow(a,3), F.ft(1,2,-1)); z4 = F.fsub(F.fpow(a,4), F.ft(-1,1,3))
print("minpoly residuals:", z3[:3], z4[:3])
x = F.ft(1,2,0,3)
print("x*inv(x):", F.fmul(x, F.finv(x))[:3])
print("norms: 1+a:", F.fnorm(F.fadd(F.ONE,a)), " 2+a:", F.fnorm(F.ft(2,1)), " a:", F.fnorm(a), " 1-a:", F.fnorm(F.fsub(F.ONE,a)), " 2:", F.fnorm(F.ft(2)))
import mpmath as mp
mp.mp.dps=40
print("feval (2-a^2)/2:", F.feval(F.ft(2,0,-1,2)), "expect 0.222520933956314404289")
# division test
y = F.fdiv(F.ft(1,1), F.ft(2,0,-1,2))
print("(1+a)/h:", y, " eval:", F.feval(y))
print("check *h:", F.fmul(y, F.ft(2,0,-1,2))[:3])
'''
import subprocess
r = subprocess.run(['<project>/.venv/bin/python','-c',test],capture_output=True,text=True)
print(r.stdout); print(r.stderr[-500:] if r.stderr else '')
