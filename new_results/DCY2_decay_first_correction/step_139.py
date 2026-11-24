from mpmath import mp, mpf, sqrt, exp
mp.dps = 50
def v(u):
    return sqrt((1-exp(-2*u))/(2*u)) if u != 0 else mpf(1)
def rho3(u):
    return v(u) - 1 + u/2 - 5*u**2/24
mx = 0; mxu = None
for i in range(1, 2001):
    u = mpf(i)/2000
    val = abs(rho3(u))/u**3
    if val > mx: mx = val; mxu = u
print("max |rho3|/u^3 on (0,1]:", mp.nstr(mx, 12), "at u =", mxu)
# also check what rho3(u)/u^3 looks like near u=0.5 and u=1
for u in [mpf('0.25'), mpf('0.5'), mpf('0.75'), mpf('1.0')]:
    print(f"  u={u}: |rho3|/u^3 = {mp.nstr(abs(rho3(u))/u**3, 8)}")
# what does rho3 look like: v = 1 - u/2 + 5u^2/24 - u^3/16 + 79u^4/5760 - ...
# rho3 = v - 1 + u/2 - 5u^2/24 = -u^3/16 + 79u^4/5760 - 3u^5/1280 + ...
# |rho3|/u^3 = |-1/16 + 79u/5760 - 3u^2/1280 + ...|
print("rho3 series leading: -u^3/16 + 79u^4/5760 - 3u^5/1280")
print("|-1/16 + 79u/5760 - 3u^2/1280| at u=0.5:", abs(-mpf(1)/16 + 79*mpf('0.5')/5760 - 3*mpf('0.5')**2/1280))
print("|-1/16 + 79u/5760 - 3u^2/1280| at u=1:", abs(-mpf(1)/16 + 79*mpf('1')/5760 - 3*mpf('1')**2/1280))
