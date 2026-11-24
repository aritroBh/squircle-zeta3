from mpmath import mp, mpf, sqrt, exp
mp.dps = 50
def v(u):
    return sqrt((1-exp(-2*u))/(2*u)) if u != 0 else mpf(1)
def rho3(u):
    return v(u) - 1 + u/2 - 5*u**2/24
# where is the max of |rho3|/u^3 on (0,1]?
mx = 0; mxu = None
for i in range(1, 4001):
    u = mpf(i)/4000
    val = abs(rho3(u))/u**3
    if val > mx: mx = val; mxu = u
print("max |rho3|/u^3 on (0,1]:", mp.nstr(mx, 12), "at u =", mxu)
# rho3 series: -u^3/16 + 79u^4/5760 - 3u^5/1280 + 71u^6/193536 - 113u^7/1935360 + 3053u^8/464486400
# rho3/u^3 = -1/16 + 79u/5760 - 3u^2/1280 + 71u^3/193536 - 113u^4/1935360 + 3053u^5/464486400
# check the sign of rho3/u^3 (alternating?)
for u in [mpf('0.1'), mpf('0.3'), mpf('0.5'), mpf('0.7'), mpf('0.9'), mpf('1.0')]:
    print(f"  u={u}: rho3/u^3 = {mp.nstr(rho3(u)/u**3, 10)}")
# tail bound: |rho3(u)/u^3 + 1/16 - 79u/5760 + 3u^2/1280| <= 71u^3/193536 + 113u^4/1935360 + tail
# the alternating structure: rho3/u^3 = -1/16 + 79u/5760 - 3u^2/1280 + ...
# with c6 = +71/193536, c7 = -113/1935360, c8 = +3053/464486400
# is rho3/u^3 between -1/16+79u/5760-3u^2/1280 and -1/16+79u/5760?
for u in [mpf('0.2'), mpf('0.5'), mpf('0.8'), mpf('1.0')]:
    lo = -mpf(1)/16 + 79*u/5760 - 3*u**2/1280
    hi = -mpf(1)/16 + 79*u/5760
    print(f"  u={u}: rho3/u^3 = {mp.nstr(rho3(u)/u**3, 8)}, in [{mp.nstr(lo,8)}, {mp.nstr(hi,8)}]? {lo <= rho3(u)/u**3 <= hi}")
