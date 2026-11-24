import mpmath as mp
mp.mp.dps = 80
r2=mp.sqrt(2); mpp=mp.pi
def psi(t): return r2*t - mp.log(1+r2*t+t*t)
def G(y, sg):
    rho = 1+y*y*sg; t = rho*y
    return psi(t) - (r2/3)*t**3/(1+r2*t+t*t) - r2*t*(rho-1)**2/rho**2
y100 = mp.sqrt(mpp/200)
print(f"y(N=100) = {mp.nstr(y100,15)} ; t_max = y(1+y^2) = {mp.nstr(y100*(1+y100**2),15)} (claim 0.12731 < 0.13)")
# dense scan over (y, sg): y in (0, y100], sg in [0,1]
mn = mp.mpf('9'); at=None
NY, NS = 400, 400
for iy in range(1, NY+1):
    y = y100*iy/NY
    for isg in range(NS+1):
        sg = mp.mpf(isg)/NS
        v = G(y, sg)
        if v < mn: mn=v; at=(y,sg)
print(f"min margin/(N/2) = G(y,sg) over full domain y in (0,y(100)], sg in [0,1] ({NY}x{NS} grid) = {mp.nstr(mn,12)}")
print(f"   at y = {mp.nstr(at[0],15)} (N = {mp.nstr(mpp/(2*at[0]**2),10)}), sigma = {mp.nstr(at[1],6)}")
# scan finer around the argmin
y0c, s0c = at
mn2 = mp.mpf('9'); at2=None
for iy in range(-50, 51):
    y = y0c*(1+mp.mpf(iy)/5000)
    for isg in range(-50, 51):
        sg = s0c + mp.mpf(isg)/5000
        if sg < 0 or sg > 1: continue
        v = G(y, sg)
        if v < mn2: mn2=v; at2=(y,sg)
print(f"refined min G = {mp.nstr(mn2,12)} at y={mp.nstr(at2[0],15)}, sigma={mp.nstr(at2[1],8)}")
