import subprocess, mpmath as mp
mp.mp.dps = 80
r2=mp.sqrt(2); mpp=mp.pi
def psi(t): return r2*t - mp.log(1+r2*t+t*t)
def G(y, sg):
    rho = 1+y*y*sg; t = rho*y
    return psi(t) - (r2/3)*t**3/(1+r2*t+t*t) - r2*t*(rho-1)**2/rho**2
y100 = mp.sqrt(mpp/200)
# (d1) G(y,1) zero crossing (R5: y=0.109634501782, N=130.68489)
root = mp.findroot(lambda yy: G(yy,1), (mp.mpf('0.109'), mp.mpf('0.111')))
print(f"(d1) tight-box failure at sigma=1: G(y,1)=0 at y={mp.nstr(root,15)}, N={mp.nstr(mpp/(2*root**2),10)}")
mn = mp.mpf('9'); at=None
for i in range(0,401):
    y = y100*(1-mp.mpf(i)/400*0.3)   # scan y in [0.7*y100, y100] where min lies
    v = G(y,1)
    if v < mn: mn=v; at=y
print(f"     min G(y,1) on y in [0.7,1]*y100 = {mp.nstr(mn,12)} at N={mp.nstr(mpp/(2*at**2),8)}  (R5: -5.4722904e-6 at y(100))")
# (d2) min G over sigma in [0,0.5] (saddle side), y in (0,y100]
mn2 = mp.mpf('9'); at2=None
for iy in range(1, 201):
    y = y100*iy/200
    for isg in range(0, 101):
        sg = mp.mpf(isg)/200
        v = G(y, sg)
        if v < mn2: mn2=v; at2=(y,sg)
print(f"(d2) min G over sigma in [0,0.5], y<=y100 = {mp.nstr(mn2,12)} at y={mp.nstr(at2[0],10)}, sig={mp.nstr(at2[1],6)}  (R5: +2.56e-14)")
# (d3) tight-box margins at saddle itself (sigma=sigma*(y)) for N in [100,400]
def saddle_sig(yv):
    yv=mp.mpf(yv); f=lambda tt: tt/mp.sqrt(1+tt**2/(1+r2*tt))-yv
    t=mp.findroot(f, yv*(1+yv**2/2)); rho=mp.sqrt(1+t**2/(1+r2*t))
    return (rho-1)/yv**2
mn3 = mp.mpf('9'); at3=None
for i in range(0, 2001):
    Nv = 100 + 300*i/2000; y = mp.sqrt(mpp/(2*Nv)); sg = saddle_sig(y)
    v = G(y, sg)
    if v < mn3: mn3=v; at3=Nv
print(f"(d3) tight-box margin G(y,sigma*(y)) on saddle, N in [100,400]: min = {mp.nstr(mn3,12)} at N={mp.nstr(at3,8)}")
print()
# (e) run certify_H3.py verbatim
r = subprocess.run(['<project>/.venv/bin/python',
    '<project>/research/08_theory/D13s_fix_recertify/D13s_work/fixed/certify_H3.py'],
    capture_output=True, text=True, cwd='/tmp')
print("(e) certify_H3.py output:"); print(r.stdout); print(r.stderr[-500:] if r.stderr else "")
