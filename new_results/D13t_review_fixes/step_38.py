import mpmath as mp
mp.mp.dps = 80
r2=mp.sqrt(2); mpp=mp.pi
def psi(t): return r2*t - mp.log(1+r2*t+t*t)
def G(y, sg):
    rho = 1+y*y*sg; t = rho*y
    return psi(t) - (r2/3)*t**3/(1+r2*t+t*t) - r2*t*(rho-1)**2/rho**2
y100 = mp.sqrt(mpp/200)
# (a) G(y,1) sign: scan
print("(a) G(y, sigma=1) at selected y (N):")
for Nv in [100, 105, 110, 120, 150, 200, 400, 1000, 1e4]:
    y = mp.sqrt(mpp/(2*mp.mpf(Nv)))
    print(f"   N={Nv:>7}: G(y,1) = {mp.nstr(G(y,1),8)}   G(y,0.5) = {mp.nstr(G(y,mp.mpf('0.5')),8)}   G(y,0.75) = {mp.nstr(G(y,mp.mpf('0.75')),8)}")
# where does G(y,1) cross zero near N=100?
f = lambda yy: G(yy, 1)
root = mp.findroot(f, (mp.mpf('0.125'), mp.mpf('0.127')))
print(f"   G(y,1)=0 at y = {mp.nstr(root,12)}  i.e. N = {mp.nstr(mpp/(2*root**2),8)}")
# (b) min of G over sigma in [0, 0.5], y in (0, y100]
mn = mp.mpf('9'); at=None
for iy in range(1, 201):
    y = y100*iy/200
    for isg in range(0, 101):
        sg = mp.mpf(isg)/200   # sigma in [0, 0.5]
        v = G(y, sg)
        if v < mn: mn=v; at=(y,sg)
print(f"(b) min G over y in (0,y100], sigma in [0,0.5] = {mp.nstr(mn,10)} at y={mp.nstr(at[0],10)} (N={mp.nstr(mpp/(2*at[0]**2),8)}), sigma={mp.nstr(at[1],6)}")
