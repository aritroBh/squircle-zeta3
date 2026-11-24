import subprocess
PY = '<project>/.venv/bin/python'
chk = r'''
import mpmath as mp
mp.mp.dps = 50
def Jx(c):
    P = -1/c; Q = -(1-c)/(2*c); Q2 = -(1-c)/c
    H1 = -mp.log(1-c)
    ell = lambda x: mp.log(x) + mp.log(1-x/P) - mp.log(1-x/Q) + H1
    emm = lambda x: mp.log(1-x) + mp.log(1-x/Q2) - mp.log(1-x/Q)
    fP = lambda x: 1/(x-P); fQ = lambda x: 1/(x-Q); fQ2 = lambda x: 1/(x-Q2)
    def g(x):
        return ell(x)**2*(fP(x)-fQ2(x)) + ell(x)*emm(x)*(fQ(x)-2*fP(x))
    return mp.quad(g, [0, mp.mpf(1)/2, 1])
def Jval(c):
    r = (1-c)/(1+c)
    def lt(w): return mp.log((w-r)*(w+1)/(2*(1-r)*w))
    def lmt(w): return mp.log((1-w)*(w+r)/(2*(1-r)*w))
    def f(w):
        return (1/(w+1) - 1/(w+r))*lt(w)**2 + (1/w - 2/(w+1))*lt(w)*lmt(w)
    return mp.quad(f, [r, mp.sqrt(r), 1])
for c in (mp.mpf(1)/2, mp.mpf(3)/5, mp.mpf(1)/3):
    a = Jx(c); b = Jval(c)
    print(f"c={mp.nstr(c,8)}: Jx={mp.nstr(a,25)} Jw={mp.nstr(b,25)} diff={mp.nstr(abs(a-b),3)}")
'''
r = subprocess.run(['nice','-n','10',PY,'-c',chk], capture_output=True, text=True)
print(r.stdout)
