import mpmath as mp
mp.mp.dps = 50
def G_int(word, z=1, dps=50):
    mp.mp.dps = dps
    if len(word) == 0: return mp.mpf(1)
    a = word[-1]; sub = word[:-1]
    pts = [mp.mpf(0)]
    if 0 < a < z: pts.append(mp.mpf(a))
    pts.append(mp.mpf(z))
    f = lambda t: G_int(sub, t, dps)/(t-a)
    return mp.quad(f, pts)
g23 = G_int((2,3))
cand = -mp.log(2)**2 - mp.re(mp.polylog(2,-2))
print("G(2,3;1):", g23)
print("cand     :", cand)
print("diff (should be c*pi^2):", (g23-cand)/(mp.pi**2/6))
