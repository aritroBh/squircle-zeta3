# SM6: verify the reduction c6, c7 in gamma_ab against quadrature, and against Richardson of item 1.
import mpmath as mp
mp.mp.dps = 80
z = mp.zeta
def gam(a, b):
    f = lambda s: s**a*mp.log(1+mp.exp(-2*s))**b
    pts = [0, 1, 2, 4, 8, 16, 32, 64, 128]
    return mp.quad(f, pts)/(mp.factorial(a)*mp.factorial(b))
g = {}
for (a, b) in [(0,1),(0,2),(1,1),(0,3),(1,2),(2,1),(0,4),(1,3),(3,1),(0,5),(1,4),(3,2),(0,6),(1,5),(3,3),(5,1),(2,2),(2,3),(2,4),(4,1),(4,2),(1,5)]:
    g[a, b] = gam(a, b)
c = {}
c[4] = 8*(g[0,3]+g[1,2]); c[5] = 8*(g[0,4]+g[1,3]+3*g[3,1])
c[6] = 8*(g[0,5]+g[1,4]+3*g[3,2]); c[7] = 8*(g[0,6]+g[1,5]+3*g[3,3]+45*g[5,1])
print("c4 =", mp.nstr(c[4], 40), " zeta4/4 =", mp.nstr(z(4)/4, 40))
print("c5 =", mp.nstr(c[5], 40), " closed:", mp.nstr(mp.mpf(107)/32*z(5)-z(2)*z(3), 40))
print("c6 (gamma formula) =", mp.nstr(c[6], 50))
print("c7 (gamma formula) =", mp.nstr(c[7], 50))
# reduced forms
c6r = 40*g[3,2] + mp.mpf(11)/16*z(6) - z(3)**2/2
c7r = 40*g[3,3] + mp.pi**4*z(3)/120 + mp.mpf(2339)/512*z(7)
print("c6 reduced  =", mp.nstr(c6r, 50), " diff", mp.nstr(c6r-c[6], 3))
print("c7 reduced  =", mp.nstr(c7r, 50), " diff", mp.nstr(c7r-c[7], 3))
c6R = mp.mpf("0.042950827880284033929433804835135265510296727692058")
print("c6 Richardson =", mp.nstr(c6R, 50), " minus formula:", mp.nstr(c6R-c[6], 3))
# verify each solved relation numerically
rel = {
 (0,3): -g[1,2] + mp.pi**4/2880,
 (0,4): -g[1,3] - mp.pi**2*z(3)/48 + 31*z(5)/128,
 (0,5): -g[1,4] + 2*g[3,2] + 11*mp.pi**6/120960 - z(3)**2/16,
 (0,6): -g[1,5] + 2*g[3,3] + mp.pi**4*z(3)/960 - 31*z(7)/256,
 (2,2): mp.pi**2*z(3)/96 - 29*z(5)/256,
 (2,3): -3*g[3,2] - mp.pi**6/24192 + z(3)**2/32,
 (2,4): -3*g[3,3] - mp.pi**4*z(3)/576 - mp.pi**2*z(5)/96 + 315*z(7)/1024,
 (3,1): 15*z(5)/256, (4,1): 31*mp.pi**6/967680, (5,1): 63*z(7)/4096,
 (4,2): 7*mp.pi**4*z(3)/23040 + mp.pi**2*z(5)/384 - 251*z(7)/4096,
}
for k, v in rel.items():
    print("relation gamma%d%d: residual %s" % (k[0], k[1], mp.nstr(g[k]-v, 3)))
print("g32 =", mp.nstr(g[3,2], 50)); print("g33 =", mp.nstr(g[3,3], 50))
