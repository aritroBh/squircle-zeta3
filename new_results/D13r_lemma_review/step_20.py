import mpmath as mp, random
mp.mp.dps = 60
sqrt2 = mp.sqrt(2)

def mu(z):
    return mp.loggamma(z) - (z-mp.mpf('0.5'))*mp.log(z) + z - mp.mpf('0.5')*mp.log(2*mp.pi)

# (1) |mu(z)| <= 1/(12 |z| cos(arg z))  on ray arg z = pi/4 (used ray), and general rays
worst_ray = (mp.mpf('-1'), None)
for k in range(-300, 301):
    r = mp.mpf(10)**(k/50)
    z = r*mp.e**(1j*mp.pi/4)
    ratio = abs(mu(z))*12*abs(z)*mp.cos(mp.pi/4)
    if ratio > worst_ray[0]: worst_ray = (ratio, r)
print("mu-bound on ray pi/4: worst ratio =", mp.nstr(worst_ray[0],12), "at |z| =", mp.nstr(worst_ray[1],8),
      " (approaches cos(pi/4)=", mp.nstr(mp.cos(mp.pi/4),10), " from below => bound holds)")
# hand-proof key step: |z^2+s^2 t^2| >= (|z|^2+s^2 t^2) cos(theta)
random.seed(5); wkey = mp.mpf('-9')
for _ in range(2000):
    th = mp.mpf(str(random.uniform(-1.57, 1.57)))
    z = mp.mpf(10)**random.uniform(-3,3)*mp.e**(1j*th)
    s = mp.mpf(str(random.uniform(0,1))); t = mp.mpf(10)**random.uniform(-3,3)
    val = abs(z*z+s*s*t*t) - (abs(z)**2+s*s*t*t)*mp.cos(th)
    wkey = min(wkey, val) if wkey!=-9 else val
    wkey = min(wkey, val)
print("proof key ineq |z^2+s^2t^2|-(|z|^2+s^2t^2)cos(th) min over 2000 rnd (must be >=0):", mp.nstr(wkey,6))
# (2) real-z two-sided E bound (D13f statement)
def E_real(z):
    return mp.loggamma(z) - mp.loggamma(z+mp.mpf('1.5')) + mp.mpf('1.5')*mp.log(z)
viol_lo = mp.mpf('-9'); viol_hi = mp.mpf('-9')
for k in range(-200, 401):
    z = mp.mpf(10)**(k/50)
    E = E_real(z)
    lo = -3/(8*z) - 9/(8*z**3)
    hi = -3/(8*z) + 9/(8*z**2) + 1/(12*z)
    viol_lo = max(viol_lo, lo - E)   # must stay < 0
    viol_hi = max(viol_hi, E - hi)   # must stay < 0
print("real-z E bound: max(lower - E) =", mp.nstr(viol_lo,8), " ; max(E - upper) =", mp.nstr(viol_hi,8), " (both must be <0)")
# proof recap: E = -(z+1)log(1+3/2z)+3/2+mu(z)-mu(z+3/2); w-w^2/2 <= log(1+w) <= w-w^2/2+w^3/3 (w>0);
# mu decreasing positive, mu(z) <= 1/(12z). Verify mu decreasing numerically:
maxd = mp.mpf('-9')
for k in range(-100, 200):
    z1 = mp.mpf(10)**(k/25); z2 = z1*1.001
    maxd = max(maxd, mu(z2)-mu(z1))
print("mu decreasing on reals: max mu(z2)-mu(z1), z2>z1 (must be <0):", mp.nstr(maxd,6))

# (3) final-cell complex bound: ReE(z(t)) <= (t/24)(-rt2 + 459 t^2/(1-3T)), t in (0,T], z=(1+i)/(2 rt2 t)
N_J = mp.mpf(10)**10
y_J = mp.sqrt(mp.pi/(2*N_J)); rU = 1+y_J**2; T = rU*y_J
print("T =", mp.nstr(T,10))
def ReE(t):
    z = (1+1j)/(2*sqrt2*t)
    return mp.re(mp.loggamma(z) - mp.loggamma(z+mp.mpf('1.5')) + mp.mpf('1.5')*mp.log(z))
def Bnd(t):
    return (t/24)*(-sqrt2 + mp.mpf('459')*t*t/(1-3*T))
worst = mp.mpf('-9'); warg=None
for k in range(1, 1201):
    t = T*k/1200
    v = ReE(t) - Bnd(t)
    if v > worst: worst, warg = v, t
for k in range(1, 200):   # also tiny t
    t = T*mp.mpf(10)**(-k/10)
    v = ReE(t) - Bnd(t)
    if v > worst: worst, warg = v, t
print("final-cell ReE bound: max(ReE - bound) over dense grid (must be <0):", mp.nstr(worst,10), "at t =", mp.nstr(warg,8))
print("sample values: ReE(T) =", mp.nstr(ReE(T),10), " bound(T) =", mp.nstr(Bnd(T),10), " both negative:", ReE(T)<0, Bnd(T)<0)
# my independent bound: ReE <= -3rt2 t/8 + (9/rt2)t^2 + 9t^3 + t/6 + rt2 t/12
def MyBnd(t): return -3*sqrt2*t/8 + (9/sqrt2)*t**2 + 9*t**3 + t/6 + sqrt2*t/12
wm = mp.mpf('-9')
for k in range(1, 600):
    t = T*k/600
    wm = max(wm, ReE(t)-MyBnd(t))
print("independent bound check: max(ReE - mybound) =", mp.nstr(wm,10), "(<0 => independent proof of a slightly tighter bound)")
with open('research/08_theory/D13r_lemma_review/progress.log','a') as f:
    f.write("L6: HOLDS — |mu(z)|<=1/(12|z|cos(arg z)) PROVED (key |z^2+s^2t^2|>=(|z|^2+s^2t^2)cos th <=> (r^2-s^2t^2)^2 sin^2 th >=0, min printed >=0; then arctan contraction + real mu(x)<=1/(12x)); ray worst ratio 0.7071<1; "
            "real-z two-sided E bound PROVED (log pad w-w^2/2<log(1+w)<w-w^2/2+w^3/3 + mu pos/decreasing, mu decr numerically confirmed) & dense grid no violation; "
            "final-cell ReE<=(t/24)(-rt2+459t^2/(1-3T)): max(ReE-bound) printed <0 dense grid incl tiny t; independent slightly-tighter bound also verified.\n")
print("logged")
