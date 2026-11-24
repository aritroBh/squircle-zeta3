import mpmath as mp
mp.mp.dps = 60
sqrt2 = mp.sqrt(2); sqpi = mp.sqrt(mp.pi); a = mp.pi/(2*sqrt2)
b = mp.e**(1j*mp.pi/4)

def Hc(s): return sqpi/4*mp.gamma(s/2)/mp.gamma(s/2+mp.mpf('1.5'))
def g(x): return Hc(b*x)/x
# gamma = log g ; gamma' and gamma'' via digamma
def gp(x):
    w = b*x/2
    return (b/2)*(mp.psi(0,w) - mp.psi(0,w+mp.mpf('1.5'))) - 1/x
def gpp(x):
    w = b*x/2
    return (b**2/4)*(mp.psi(1,w) - mp.psi(1,w+mp.mpf('1.5'))) + 1/x**2
def G1b(x):
    t = x/(2*sqrt2)
    return mp.mpf('0.75')*(1/t**2 + 1/t) + 1/x
def G2b(x):
    t = x/(2*sqrt2)
    return mp.mpf('0.375')*(2/t**3 + 1/t**2) + 1/x**2

# cross-check gp/gpp against finite differences of log g
x0 = mp.mpf('7.9')
h = mp.mpf('1e-12')
fd1 = (mp.log(g(x0+h))-mp.log(g(x0-h)))/(2*h)
fd2 = (mp.log(g(x0+h))-2*mp.log(g(x0))+mp.log(g(x0-h)))/h**2
print("def check: gp vs fd:", abs(gp(x0)-fd1), " gpp vs fd:", abs(gpp(x0)-fd2))

worst1 = (mp.mpf('2'),None); worst2 = (mp.mpf('2'),None)
import random
random.seed(13)
xs = []
for k in range(-60, 200):
    xs.append(mp.mpf(10)**(k/10))
xs += [mp.mpf(str(random.uniform(0.01, 1e7))) for _ in range(300)]
xs += [mp.mpf('5.3'), mp.mpf('7.98'), mp.mpf('40'), mp.mpf('400'), mp.mpf('4e4')]
for x in xs:
    r1 = abs(gp(x))/G1b(x); r2 = abs(gpp(x))/G2b(x)
    if r1 > worst1[0]: worst1 = (r1, x)
    if r2 > worst2[0]: worst2 = (r2, x)
print("L1 gamma': worst ratio |g'|/bound =", mp.nstr(worst1[0],8), "at x =", mp.nstr(worst1[1],8))
print("L1 gamma'': worst ratio |g''|/bound =", mp.nstr(worst2[0],8), "at x =", mp.nstr(worst2[1],8))

# auxiliary inequalities: psi'(t) <= 1/t^2+1/t ; |psi''(t)| <= 2/t^3+1/t^2  (t real >0)
w3=(mp.mpf('2'),None); w4=(mp.mpf('2'),None)
for k in range(-80, 120):
    t = mp.mpf(10)**(k/8)
    r3 = mp.psi(1,t)/(1/t**2+1/t)
    r4 = abs(mp.psi(2,t))/(2/t**3+1/t**2)
    if r3>w3[0]: w3=(r3,t)
    if r4>w4[0]: w4=(r4,t)
print("psi' bound worst ratio:", mp.nstr(w3[0],8), "at t=", mp.nstr(w3[1],6))
print("psi'' bound worst ratio:", mp.nstr(w4[0],8), "at t=", mp.nstr(w4[1],6))

# series identity: psi(w)-psi(w+3/2) = -(3/2) sum 1/((n+w)(n+w+3/2)), w=t(1+i)
t = mp.mpf('1.2345'); w = t*(1+1j)
s = mp.nsum(lambda n: 1/((n+w)*(n+w+mp.mpf('1.5'))), [0, mp.inf])
res = (mp.psi(0,w)-mp.psi(0,w+mp.mpf('1.5'))) + mp.mpf('1.5')*s
print("series identity residual:", abs(res))
# where are the maxima of the ratios? refine near worst points
for name, f, B in [("g'",gp,G1b),("g''",gpp,G2b)]:
    best=(mp.mpf('0'),None)
    for k in range(-400,400):
        x = mp.mpf(10)**(k/100)
        r = abs(f(x))/B(x)
        if r>best[0]: best=(r,x)
    print(f"L1 {name}: dense-grid worst ratio = {mp.nstr(best[0],10)} at x={mp.nstr(best[1],10)}")
