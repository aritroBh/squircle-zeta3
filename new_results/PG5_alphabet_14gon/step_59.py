chk2 = r'''
import mpmath as mp, pickle, sys
sys.path.insert(0,'.')
import ftfield as F
mp.mp.dps = 60
m = 7
g = pickle.load(open('geom14.pkl','rb'))
E3 = pickle.load(open('E3n14.pkl','rb'))['E3n']
data = g['data']; ORDER = [3,4,1,2,5,6,0]
apn = [[F.feval(data[i][0],60), F.feval(data[i][1],60)] for i in ORDER]
hh = F.feval(g['h'],60); LL = mp.sqrt(F.feval(g['L'],60))
def E3eval(t):
    logs = [mp.log(apn[i][0]+apn[i][1]*t) for i in range(7)]
    val = mp.mpf(0)
    for k,v in E3.items():
        poly = ((F.feval(v[3],60)*t+F.feval(v[2],60))*t+F.feval(v[1],60))*t+F.feval(v[0],60)
        lv = mp.mpf(1)
        for i,e in enumerate(k):
            if e: lv *= logs[i]**e
        val += poly*lv
    return val
def r_and_rp(th, p):
    f = sum(mp.exp(p*mp.log(abs(mp.cos(th - i*mp.pi/m)))) for i in range(m))
    fp = sum(-mp.sin(th - i*mp.pi/m)*mp.sign(mp.cos(th-i*mp.pi/m))*p*mp.exp((p-1)*mp.log(abs(mp.cos(th-i*mp.pi/m)))) for i in range(m))
    r = f**(-1/p); rp = -(1/p)*f**(-1/p-1)*fp
    return r, rp
def Jfun(t, d):
    y = (t-mp.mpf('0.5'))*LL
    th = mp.atan(y/hh)
    thp = (LL/hh)/(1+(y/hh)**2)
    r, rp = r_and_rp(th, 1+d)
    return mp.sqrt(r*r+rp*rp)*thp/LL
def d3J(t, N=12, rho=mp.mpf('0.03')):
    tot = mp.mpc(0)
    for j in range(N):
        z = mp.exp(2*mp.pi*1j*j/N)
        tot += z**(-3)*Jfun(t, rho*z)
    return 6*tot/N/rho**3
# integral from the independent side: 14*L*int d3J/6
Ij, err = mp.quad(lambda t: d3J(t)/6, [0,1], error=True)
print("int d3J/6 =", mp.nstr(Ij,25), " err", mp.nstr(err,3))
gold = pickle.load(open('gold14.pkl','rb'))
P3g = mp.mpf(gold['P3_gold'])
print("14*L*that =", mp.nstr(14*LL*Ij,25))
print("gold      =", mp.nstr(P3g,25))
print("diff:", mp.nstr(abs(14*LL*Ij-P3g),4))
# now E3eval integral by Simpson on fine grid (avoid endpoint issues)
import numpy as np
def simpson(f, NN=20000):
    h = mp.mpf(1)/NN
    s = f(mp.mpf('0.5')*h)
    tot = f(h/2)*4
    for k in range(1,NN):
        x = k*h
        tot += (2 if k%2==0 else 4)*f(x) if False else (4 if k%2 else 2)*f(x)
    tot += f(1-h/2)*4 + f(1-h/2)
    return None
# simpler: tanh-sinh with more precision on subintervals
I1 = mp.quad(E3eval, [0, 0.5])
I2 = mp.quad(E3eval, [0.5, 1])
print("int E3eval split:", mp.nstr(I1+I2,25))
print("int d3J/6       :", mp.nstr(Ij,25))
print("diff of integrals:", mp.nstr(abs(I1+I2-Ij),4))
# midpoint rule crude
N=4000
s = sum(E3eval((k+mp.mpf('0.5'))/N) for k in range(N))/N
print("midpoint int E3eval:", mp.nstr(s,20))
'''
open('step_08b_intcheck.py','w').write(chk2)
import subprocess, time
t0=time.time()
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_08b_intcheck.py'],capture_output=True,text=True,timeout=1500)
print("elapsed", time.time()-t0)
print(r.stdout); print(r.stderr[-800:] if r.stderr else '')
