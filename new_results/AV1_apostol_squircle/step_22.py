code = r'''
# step_19.py : d1 to 85 digits via smooth transformed integrals; redo d9 residual-Richardson.
import mpmath as mp
mp.mp.dps = 130
delta = mp.mpf('0.5')
def L(u): return mp.log(1-mp.exp(-u))
# near piece [0,delta]: u=e^{-v}, v in [log(1/delta), inf)
def near_integrand(v):
    u = mp.exp(-v)
    w = u - L(u)          # = -Phi > 0 ; as v->inf, w ~ v + 1.5 e^{-v}
    return (mp.log(w) + v)*u
near = mp.quad(near_integrand, [mp.log(1/delta), 1, 3, 8, 20, 50, mp.inf])
far = mp.quad(lambda u: mp.log((u-L(u))/u), [delta, 1, 2, 5, 15, 40, mp.inf])
d1 = near+far
print('d1 =', mp.nstr(d1, 60))

# check via Richardson on series
def A_series(p, tol=mp.mpf('1e-125')):
    p = mp.mpf(p); s = mp.mpf(0); m = 0
    while True:
        m += 1
        t = mp.gamma(1+m/p)**2/(m**2*mp.gamma(1+2*m/p))
        s += t
        if m > 30 and abs(t) < tol: break
    return s
ps = list(range(80, 801, 20))
xs=[]; ys=[]
for p in ps:
    xs.append(mp.mpf(1)/p); ys.append((mp.zeta(2)-A_series(p))*p)
deg=9
n=len(xs); M=mp.matrix(n,deg)
for i,x in enumerate(xs):
    for j in range(deg): M[i,j]=x**j
b=mp.matrix(n,1)
for i,y in enumerate(ys): b[i,0]=y
sol=mp.lu_solve(M.T*M,M.T*b)
print('Richardson d1:', mp.nstr(sol[0,0],30), ' |diff|', mp.nstr(abs(sol[0,0]-d1),3))

# redo d9 residual Richardson with corrected d1
mp.mp.dps = 105
d2 = -mp.zeta(2)/2; d3 = mp.zeta(3)/6
d5 = (mp.zeta(2)*mp.zeta(3)-3*mp.zeta(5))/60
d7 = mp.zeta(7)/14 - mp.zeta(2)*mp.zeta(5)/42 - mp.zeta(2)**2*mp.zeta(3)/140
d9cand = (-428400*mp.zeta(9)+136080*mp.zeta(2)*mp.zeta(7)+79*mp.pi**6*mp.zeta(3)+1134*mp.pi**4*mp.zeta(5)-10080*mp.zeta(3)**3)/1814400
xs=[]; ys=[]
for p in list(range(60, 601, 30)):
    e = mp.mpf(1)/p
    D = mp.zeta(2) - A_series(p, tol=mp.mpf('1e-100'))
    R = (D - d1*e - d2*e**2 - d3*e**3 - d5*e**5 - d7*e**7)/e**9
    xs.append(e); ys.append(R)
deg=4
n=len(xs); M=mp.matrix(n,deg)
for i,x in enumerate(xs):
    for j in range(deg): M[i,j]=x**(2*j)
b=mp.matrix(n,1)
for i,y in enumerate(ys): b[i,0]=y
sol=mp.lu_solve(M.T*M,M.T*b)
print('d9 residual-Richardson:', mp.nstr(sol[0,0],15), ' candidate:', mp.nstr(d9cand,15), ' |diff|', mp.nstr(abs(sol[0,0]-d9cand),3))
print('d1 full:', mp.nstr(d1, 90))
'''
open('steps/step_19.py','w').write(code)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','steps/step_19.py'],capture_output=True,text=True,timeout=2400)
print(r.stdout); print(r.stderr[-800:] if r.returncode else '')
open('progress.log','a').write('step_19: d1 recomputed to 85+ digits (old 45-digit value wrong past 1e-11)\n')
