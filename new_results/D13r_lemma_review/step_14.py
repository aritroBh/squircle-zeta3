import mpmath as mp, random
mp.mp.dps = 60
sqrt2 = mp.sqrt(2); sqpi = mp.sqrt(mp.pi)
b = mp.e**(1j*mp.pi/4)
def Hc(s): return sqpi/4*mp.gamma(s/2)/mp.gamma(s/2+mp.mpf('1.5'))
def g(x): return Hc(b*x)/x
def gp(x):
    w = b*x/2
    return (b/2)*(mp.psi(0,w) - mp.psi(0,w+mp.mpf('1.5'))) - 1/x
def gpp(x):
    w = b*x/2
    return (b**2/4)*(mp.psi(1,w) - mp.psi(1,w+mp.mpf('1.5'))) + 1/x**2
def G1b(x):
    t = x/(2*sqrt2); return mp.mpf('0.75')*(1/t**2 + 1/t) + 1/x
def G2b(x):
    t = x/(2*sqrt2); return mp.mpf('0.375')*(2/t**3 + 1/t**2) + 1/x**2

random.seed(99)
worst = {'g1':(mp.mpf('-1'),None), 'g2':(mp.mpf('-1'),None),
         'p1':(mp.mpf('-1'),None), 'p2':(mp.mpf('-1'),None)}
N=0
def feed(x):
    global N; N+=1
    r1 = abs(gp(x))/G1b(x); r2 = abs(gpp(x))/G2b(x)
    t = x/(2*sqrt2)
    r3 = mp.psi(1,t)/(1/t**2+1/t); r4 = abs(mp.psi(2,t))/(2/t**3+1/t**2)
    if r1>worst['g1'][0]: worst['g1']=(r1,x)
    if r2>worst['g2'][0]: worst['g2']=(r2,x)
    if r3>worst['p1'][0]: worst['p1']=(r3,t)
    if r4>worst['p2'][0]: worst['p2']=(r4,t)
# dense log grid over the full claim range (0, inf), finer in used range
for k in range(-300, 700): feed(mp.mpf(10)**(k/100))
# uniform random in the range the code uses: x in [4, 5e4] plus huge-x
for _ in range(400): feed(mp.mpf(str(random.uniform(4, 5e4))))
for _ in range(200): feed(mp.mpf(str(10**random.uniform(4, 8))))
# boundary regions x -> 0 and x tiny
for k in range(1, 60): feed(mp.mpf(10)**(-k/5))
print(f"sampled {N} points at 60 digits")
for k,(r,x) in worst.items():
    print(f"{k}: worst ratio = {mp.nstr(r,12)}  at {mp.nstr(x,10)}   (margin 1/ratio = {mp.nstr(1/r,6)})")
# analytic limits of the ratios as x->inf (explains worst values)
print("limit |g'|/G1b as x->inf:  (5/2)/(3*sqrt2/2+1) =", mp.nstr((mp.mpf('2.5'))/(3*sqrt2/2+1),10))
print("limit |g''|/G2b as x->inf: (27/8+1)/((3/8)*24+1) =", mp.nstr((mp.mpf(27)/8+1)/(mp.mpf('0.375')*24+1),10))
