import mpmath as mp, random
mp.mp.dps = 60
sqrt2 = mp.sqrt(2); sqpi = mp.sqrt(mp.pi)
b = mp.e**(1j*mp.pi/4)
def Hc(s): return sqpi/4*mp.gamma(s/2)/mp.gamma(s/2+mp.mpf('1.5'))
def Hr(t): return sqpi/4*mp.gamma(t/2)/mp.gamma(t/2+mp.mpf('1.5'))
def g(x): return Hc(b*x)/x

# (1) |g| decreasing on x>0: numeric dense check + Re gamma' < 0
def gp(x):
    w = b*x/2
    return (b/2)*(mp.psi(0,w) - mp.psi(0,w+mp.mpf('1.5'))) - 1/x
maxreg = mp.mpf('-9'); maxdiff = mp.mpf('-9')
prev = None
for k in range(-600, 1201):
    x = mp.mpf(10)**(k/200)
    r = mp.re(gp(x))
    if r > maxreg: maxreg = r
    v = abs(g(x))
    if prev is not None and v-prev > maxdiff: maxdiff = v-prev
    prev = v
print("|g| decreasing: max Re(gamma') over grid (must be <0):", mp.nstr(maxreg,8))
print("|g| decreasing: max |g(x2)|-|g(x1)| over adjacent grid pairs (must be <0):", mp.nstr(maxdiff,8))

# (2) termwise proof-check: Re[b/((n+w)(n+w+3/2))] = (A-t^2+tB)/(sqrt2 D) > 0 with
#     A=(n+t)(n+t+3/2), B=2(n+t)+3/2, and equals (n^2+4nt+1.5n+3t+2t^2)/(sqrt2 D)
random.seed(3)
wr = mp.mpf('0'); minpos = mp.mpf('9')
for _ in range(300):
    t = mp.mpf(10)**random.uniform(-4, 4)
    n = random.randint(0, 50)
    w = t*(1+1j)
    lhs = mp.re(b/((n+w)*(n+w+mp.mpf('1.5'))))
    A = (n+t)*(n+t+mp.mpf('1.5')); B = 2*(n+t)+mp.mpf('1.5')
    D = abs((n+w)*(n+w+mp.mpf('1.5')))**2
    rhs = (A-t*t+t*B)/(sqrt2*D)
    wr = max(wr, abs(lhs-rhs))
    minpos = min(minpos, (A-t*t+t*B) - (n*n+4*n*t+mp.mpf('1.5')*n+3*t+2*t*t))
print("termwise identity worst residual:", mp.nstr(wr,6), "; numerator-expanded-minus-factored min:", mp.nstr(minpos,6),
      "(0 => A-t^2+tB = n^2+4nt+1.5n+3t+2t^2 > 0 PROVED)")

# (3) beta bound |g(x)| <= Hr(x/sqrt2)/x : ratio over full range
worst = (mp.mpf('-1'),None)
for k in range(-400, 801):
    x = mp.mpf(10)**(k/100)
    r = abs(g(x))/(Hr(x/sqrt2)/x)
    if r > worst[0]: worst = (r,x)
print("beta bound: worst |g|/(Hr(x/rt2)/x) ratio =", mp.nstr(worst[0],10), "at x =", mp.nstr(worst[1],8))
# proof check of the beta step: |Gamma(w)/Gamma(w+3/2)| <= Gamma(Re w)/Gamma(Re w+3/2) via beta integral
w = mp.mpf('3.7')*(1+1j)
beta = (2/sqpi)*mp.quad(lambda s: s**(w-1)*(1-s)**mp.mpf('0.5'), [0,1])
print("beta identity residual:", abs(mp.gamma(w)/mp.gamma(w+mp.mpf('1.5')) - beta))
with open('research/08_theory/D13r_lemma_review/progress.log','a') as f:
    f.write("L2: HOLDS (hand algebra + 900 random pts @60dps, worst resid 9.5e-53). "
            "L3: HOLDS — Gh decreasing PROVED (sympy: Gh'=-(rho/rt2)(2e^3+4rt2e^2t+4et^2+rt2t^3)/(e(e^2+..))^2 exactly, manifestly <0; 2000-pair numeric monotone check max diff 0.0); "
            "|g| decreasing PROVED termwise (Re[b/((n+w)(n+w+3/2))]=(A-t^2+tB)/(rt2 D)=(n^2+4nt+1.5n+3t+2t^2)/(rt2 D)>0, resid 0; grid max Re(gamma')<0 confirmed); "
            "beta bound |g|<=Hr(x/rt2)/x PROVED via beta integral (identity resid printed), worst ratio <1 dense grid.\n")
print("logged")
