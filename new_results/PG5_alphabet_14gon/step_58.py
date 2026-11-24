chk = r'''
"""PG5 step_08: validate E3n for the 14-gon.
(a) pointwise: 6*E3n(t) vs independent d^3-derivative of the polar perimeter integrand
    J(t,d) = sqrt(r^2+r'^2)(theta(t),1+d) theta'(t)/L   (roots-of-unity filter, complex d).
(b) global: 84*L*int_0^1 E3n dt vs gold P3(1) (step_06b, two radii agree 7e-18).
"""
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
# independent polar integrand
def r_and_rp(th, p):
    f = sum(mp.exp(p*mp.log(abs(mp.cos(th - i*mp.pi/m)))) for i in range(m))
    fp = sum(-mp.sin(th - i*mp.pi/m)*mp.sign(mp.cos(th-i*mp.pi/m))*p*mp.exp((p-1)*mp.log(abs(mp.cos(th-i*mp.pi/m)))) for i in range(m))
    r = f**(-1/p)
    rp = -(1/p)*f**(-1/p-1)*fp
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
print("pointwise check 6*E3n(t) vs d^3 J:")
for ts in ['0.13','0.37','0.5','0.72','0.91']:
    t = mp.mpf(ts)
    a = 6*E3eval(t); b = d3J(t)
    print(f" t={ts}: {mp.nstr(a,25)}  vs  {mp.nstr(b,25)}  diff {mp.nstr(abs(a-b),3)}", flush=True)
# global check
I, err = mp.quad(E3eval, [0,1], error=True)
gold = pickle.load(open('gold14.pkl','rb'))
P3g = mp.mpf(gold['P3_gold'].split(' ')[0].strip('()')) if '(' in gold['P3_gold'] else mp.mpf(gold['P3_gold'])
val = 84*LL*I
print("int E3n =", mp.nstr(I,30), " quad err", mp.nstr(err,3))
print("84*L*int =", mp.nstr(val,30))
print("gold P3  =", mp.nstr(P3g,30))
print("|diff| =", mp.nstr(abs(val-P3g),5))
pickle.dump({'I':str(I),'val':str(val),'gold':str(P3g),'diff':str(abs(val-P3g))}, open('E3n14_check.pkl','wb'))
open('progress.log','a').write('PG5 step_08: E3n validated: pointwise diffs, 84L int vs gold diff %s\n'%mp.nstr(abs(val-P3g),3))
'''
open('step_08_e3check.py','w').write(chk)
import subprocess, time
t0=time.time()
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','step_08_e3check.py'],capture_output=True,text=True,timeout=1200)
print("elapsed", time.time()-t0)
print(r.stdout); print(r.stderr[-800:] if r.stderr else '')
