import mpmath as mp, random
mp.mp.dps = 60
sqrt2 = mp.sqrt(2); sqpi = mp.sqrt(mp.pi); a = mp.pi/(2*sqrt2); pi = mp.pi
bb = mp.e**(1j*mp.pi/4)
random.seed(23)

def Phi(x, N): return (-1+1j)*a*x + N*mp.log(bb*x/(1+bb*x))
def Phi2(x, N): return -N*(1+2*bb*x)/(x**2*(1+bb*x)**2)
def numR(t): return t**3 + 2*sqrt2*t**2 + 4*t + sqrt2
def curv2(y, rho, eta):
    num = sqrt2*eta**3*rho + 4*eta**2*rho**2*y + 2*sqrt2*eta*rho**3*y**2 + rho**4*y**3
    den = eta**2*(eta**2+sqrt2*eta*rho*y+rho**2*y**2)**2
    return (pi*rho**2*y/4)*num/den
def Ph3reg(y, rho, eta): return pi*rho**4*(rho**2*y**2+3*eta*rho*y+3*eta**2)/(eta**3*(eta**2+sqrt2*eta*rho*y+rho**2*y**2)**mp.mpf('1.5'))
def Ph4reg(y, rho, eta): return 3*pi*rho**5*(rho**3*y**3+4*eta*rho**2*y**2+6*eta**2*rho*y+4*eta**3)/(eta**4*(eta**2+sqrt2*eta*rho*y+rho**2*y**2)**2)
def G1reg(y, rho, eta): return 6*rho**2*y/eta**2 + (mp.mpf('1.5')*sqrt2+1)*rho/eta
def G2reg(y, rho, eta): return 12*sqrt2*rho**3*y/eta**3 + 4*rho**2/eta**2
def Phi3b(x,N): return 2*N*(1+3*x+3*x**2)/(x**3*(1+sqrt2*x+x*x)**mp.mpf('1.5'))
def Phi4b(x,N): return 6*N*(1+4*x+6*x**2+4*x**3)/(x**4*(1+sqrt2*x+x*x)**2)
def G1b(x):
    t = x/(2*sqrt2); return mp.mpf('0.75')*(1/t**2+1/t)+1/x
def G2b(x):
    t = x/(2*sqrt2); return mp.mpf('0.375')*(2/t**3+1/t**2)+1/x**2

R = {'curv2':0, 'ph3':0, 'ph4':0, 'g1':0, 'g2':0, 'c2':0, 'om':0, 'ms':0, 'h3lo':mp.mpf('-9'), 'h3hi':mp.mpf('-9'), 'tail':0, 'z':0}
def mx(k,v):
    global R
    R[k] = max(R[k], v)
for _ in range(400):
    y = mp.mpf(10)**random.uniform(-8,-0.4)
    N = pi/(2*y*y); th = pi/(sqrt2*y)
    rho = 1+y*y*random.random()
    eta = mp.mpf(10)**random.uniform(-0.5, 0.5)
    x0 = 1/(rho*y); x = eta*x0
    t = rho*y
    mx('curv2', abs(abs(mp.re(Phi2(x,N)))/2 - curv2(y,rho,eta))/curv2(y,rho,eta))
    mx('ph3', abs(Phi3b(x,N) - y**2*Ph3reg(y,rho,eta))/Phi3b(x,N))
    mx('ph4', abs(Phi4b(x,N) - y**3*Ph4reg(y,rho,eta))/Phi4b(x,N))
    mx('g1', abs(G1b(x) - y*G1reg(y,rho,eta))/G1b(x))
    mx('g2', abs(G2b(x) - y**2*G2reg(y,rho,eta))/G2b(x))
    # c2 = y c2reg (eta=1 case of curv2)
    c2 = abs(mp.re(Phi2(x0,N)))/2
    c2reg = (pi*rho**3/4)*numR(t)/(1+sqrt2*t+t**2)**2
    mx('c2', abs(c2 - y*c2reg)/c2)
    # om = y omreg
    Omn = sqrt2*rho + rho**2*y - 2*y*( (rho-1)/y**2 ) - y**3*((rho-1)/y**2)**2  # sig = (rho-1)/y^2
    omreg = (pi/(2*sqrt2))*Omn/(1+sqrt2*t+t**2)
    om = mp.im((-1+1j)*a + N*(1/x0 - bb/(1+bb*x0)))
    mx('om', abs(om - y*omreg)/abs(om))
for k in ['curv2','ph3','ph4','g1','g2','c2','om']:
    print(f"L7 {k}: worst REL residual over 400 pts @60dps = {mp.nstr(R[k],6)}")

# mS identity needs the TRUE saddle: solve y0 from y via the saddle equation, check RePhi(x0)+th = -(N/2)S(rho y)
def true_rho(y):
    f = lambda r: r**2*(1+sqrt2*r*y)/(1+sqrt2*r*y+r**2*y**2) - 1
    return mp.findroot(f, 1+y*y/2)
worst_ms = mp.mpf('0'); warg=None
for _ in range(120):
    y = mp.mpf(10)**random.uniform(-8,-0.4)
    N = pi/(2*y*y); th = pi/(sqrt2*y)
    rho = true_rho(y); x0 = 1/(rho*y); t = rho*y
    S = mp.log(1+sqrt2*t+t*t) - sqrt2*t
    lhs = mp.re(Phi(x0,N)) + th
    rhs = -(N/2)*S
    v = abs(lhs-rhs)/abs(lhs)
    if v > worst_ms: worst_ms, warg = v, y
print("L7 mS (true saddle): worst REL residual =", mp.nstr(worst_ms,6), "at y =", mp.nstr(warg,6))

# h3 two-sided: -S(t)/t^3 in (rt2/3)[1/(1+rt2 T+T^2), 1] for t in (0,T]
viol = mp.mpf('-9'); viol2 = mp.mpf('-9')
for T in [mp.mpf('1e-5'), mp.mpf('0.01'), mp.mpf('0.125'), mp.mpf('0.78'), mp.mpf('2.0')]:
    for k in range(1,400):
        t = T*k/400
        h3 = -(mp.log(1+sqrt2*t+t*t)-sqrt2*t)/t**3
        viol = max(viol, h3 - sqrt2/3)                    # must be <0
        viol2 = max(viol2, (sqrt2/3)/(1+sqrt2*T+T*T) - h3) # must be <0
print("L7 h3: max(h3 - rt2/3) =", mp.nstr(viol,8), " ; max(lower - h3) =", mp.nstr(viol2,8), " (both <0 required; tested T up to 2.0)")

# tail identity
wt = mp.mpf('0')
for _ in range(200):
    y = mp.mpf(10)**random.uniform(-8,-0.4); N = pi/(2*y*y)
    c_w = mp.mpf('1.0'); c2 = mp.mpf(10)**random.uniform(-1,2)
    w = mp.sqrt(mp.log(N)+c_w); d = w/mp.sqrt(2*c2)
    wt = max(wt, abs(mp.e**(-c2*d*d) - y*mp.sqrt(2/pi)*mp.e**(-c_w/2))/ (y*mp.sqrt(2/pi)*mp.e**(-c_w/2)))
print("L7 tail identity e^{-c2 d^2} = y sqrt(2/pi) e^{-c_w/2}: worst REL residual =", mp.nstr(wt,6))
