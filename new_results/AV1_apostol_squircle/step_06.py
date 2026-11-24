code = r'''
# step_05.py : exact I_{a,b} = int_0^inf u^a log^b(1-e^-u) du at 90 digits
# (analytic [0,delta] via g(u)=log((1-e^-u)/u) series + quad [delta,inf)),
# then d5, d7; PSLQ vs MZV bases; PSLQ of d1 vs broad pool.
import mpmath as mp
mp.mp.dps = 105

# series of g(u) = log((1-e^-u)/u) = sum c_k u^k  (compute by taylor)
def g_series(K):
    # (1-e^-u)/u = sum_{j>=0} (-u)^j/(j+1)!
    import sympy as sp
    u = sp.symbols('u')
    ex = sp.log(sum((-u)**j/sp.factorial(j+1) for j in range(K+3)))
    return [sp.Rational(sp.expand(ex).series(u,0,K+1).removeO().coeff(u,k)) for k in range(K+1)]
import sympy as sp
Kmax = 60
cs = g_series(Kmax)
c = [mp.mpf(sp.nsimplify(x)) for x in cs]

def gpow(j, kmax):
    # coefficients of g(u)^j up to u^kmax
    from functools import reduce
    poly = [mp.mpf(1)]
    for _ in range(j):
        new = [mp.mpf(0)]*(len(poly)+Kmax)
        for i,a in enumerate(poly):
            for k2 in range(1, Kmax+1):
                new[i+k2] += a*c[k2]
        poly = new
    return poly

def int_u_log(s, m, delta):
    # int_0^delta u^s (log u)^m du  = delta^{s+1} sum_{i=0}^m (-1)^{m-i} m!/i! (log delta)^i/(s+1)^{m-i+1}
    ld = mp.log(delta); tot = mp.mpf(0)
    for i in range(m+1):
        tot += (-1)**(m-i) * mp.factorial(m)/mp.factorial(i) * ld**i/(s+1)**(m-i+1)
    return delta**(s+1)*tot

def I_ab(a, b, delta=mp.mpf('0.5')):
    # int_0^inf u^a L^b, L = log u + g(u) on [0,delta]
    near = mp.mpf(0)
    gpow_cache = {}
    for j in range(0, b+1):
        pj = gpow(j, Kmax) if j>0 else [mp.mpf(1)]
        for k2, coeff in enumerate(pj):
            if coeff == 0: continue
            near += mp.binomial(b,j) * coeff * int_u_log(a+k2, b-j, delta)
    far = mp.quad(lambda u: u**a * mp.log(1-mp.exp(-u))**b, [delta, 1, 3, 8, 20, mp.inf])
    return near + far

# checks: I01=-zeta2, I11=-zeta3, I02=2zeta3, I31=-6zeta5
for (a,b,exact) in [(0,1,-mp.zeta(2)),(1,1,-mp.zeta(3)),(0,2,2*mp.zeta(3)),(3,1,-6*mp.zeta(5)),(2,1,-2*mp.zeta(4))]:
    v = I_ab(a,b)
    print(f'I_{a},{b} err vs exact: {mp.nstr(abs(v-exact),3)}')

I = {(a,b): I_ab(a,b) for (a,b) in [(3,1),(2,2),(1,3),(0,4),(5,1),(4,2),(3,3),(2,4),(1,5),(0,6)]}
d5 = -(-4*I[(3,1)] + 6*I[(2,2)] - 4*I[(1,3)] + I[(0,4)])/2880
d7 = (-6*I[(5,1)] + 15*I[(4,2)] - 20*I[(3,3)] + 15*I[(2,4)] - 6*I[(1,5)] + I[(0,6)])/181440
print('d5 =', mp.nstr(d5, 45))
print('d7 =', mp.nstr(d7, 45))

def pslq(val, basis, names, maxcoeff=10**9, tol=mp.mpf('1e-40')):
    rel = mp.pslq([val]+basis, tol=tol, maxcoeff=maxcoeff, maxsteps=10000)
    if rel:
        print('PSLQ hit:', rel, ' => ', names, 'residual ok')
    else:
        print('PSLQ no relation (maxcoeff=%g)'%maxcoeff)
    return rel

print('--- d5 vs {zeta5, zeta2*zeta3}')
pslq(d5, [mp.zeta(5), mp.zeta(2)*mp.zeta(3)], ['1','z5','z2z3'])
print('--- d7 vs {zeta8, zeta3*zeta5, zeta2*zeta3^2}')
pslq(d7, [mp.zeta(8), mp.zeta(3)*mp.zeta(5), mp.zeta(2)*mp.zeta(3)**2], ['1','z8','z3z5','z2z3^2'])

# zeta(3,5) via nsum for the 4-dim test
z35 = mp.nsum(lambda m: (mp.zeta(3) - mp.nsum(lambda k: k**-mp.mpf(3), [1,m]))/m**5, [1, mp.inf], tol=mp.mpf('1e-45'))
print('zeta(3,5) =', mp.nstr(z35,30))
print('--- d7 vs {zeta8, zeta3*zeta5, zeta2*zeta3^2, zeta(3,5)}')
pslq(d7, [mp.zeta(8), mp.zeta(3)*mp.zeta(5), mp.zeta(2)*mp.zeta(3)**2, z35], ['1','z8','z3z5','z2z3^2','z35'])

# d1 broad pool
d1 = mp.mpf('1.73300259017088785677509047603076917317936436')
print('--- d1 vs pool')
pslq(d1, [mp.pi, mp.euler, mp.log(2), mp.log(mp.pi), mp.zeta(2), mp.zeta(3)], ['1','pi','gam','log2','logpi','z2','z3'], maxcoeff=10**7)
pslq(d1, [mp.log(2)**2, mp.log(2)*mp.zeta(2), mp.zeta(2), mp.zeta(3)], ['1','log2^2','log2*z2','z2','z3'], maxcoeff=10**7)
pslq(d1, [mp.euler*mp.log(2), mp.log(2*mp.pi), mp.zeta(2), mp.log(2)], ['1','gam*log2','log2pi','z2','log2'], maxcoeff=10**7)
'''
open('steps/step_05.py','w').write(code)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','steps/step_05.py'],capture_output=True,text=True,timeout=2400)
print(r.stdout); print(r.stderr[-1500:] if r.returncode else '')
open('progress.log','a').write('step_05: I_ab exact at 90d; d5=zeta2*zeta3/2-3zeta5? PSLQ; d1 pool PSLQ\n')
