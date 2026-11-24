code = r'''
# step_16.py : TASK 5 square-end expansions of the variants, Richardson verification.
import mpmath as mp
mp.mp.dps = 105

def A_series_s(p, s, tol=None):
    p = mp.mpf(p); s = mp.mpf(s); tol = tol or mp.mpf('1e-100')
    tot = mp.mpf(0); m = 0
    while True:
        m += 1
        t = mp.gamma((m+s)/p)**2 / (p**2 * mp.gamma(1 + 2*(m+s)/p))
        tot += t
        if m > 30 and abs(t) < tol: break
    return tot

def B_series(p, tol=None):
    p = mp.mpf(p); tol = tol or mp.mpf('1e-100')
    tot = mp.mpf(0); n = -1
    while True:
        n += 1
        t = mp.gamma((2*n+1)/p)**2 / (p**2 * mp.gamma(1 + 2*(2*n+1)/p))
        tot += t
        if n > 30 and abs(t) < tol: break
    return tot

def richardson(vals_eps, deg):
    # vals_eps: list of (eps, deficit); fit deficit = sum_{j=1..deg} d_j eps^j
    xs = [e for e,_ in vals_eps]; ys = [v for _,v in vals_eps]
    n = len(xs); M = mp.matrix(n, deg)
    for i,x in enumerate(xs):
        for j in range(deg): M[i,j] = x**(j+1)
    b = mp.matrix(n,1)
    for i,y in enumerate(ys): b[i,0] = y
    sol = mp.lu_solve(M.T*M, M.T*b)
    return [sol[j,0] for j in range(deg)]

ps = list(range(60, 501, 20))
d1val = mp.mpf('1.73300259017088785677509047603076917317936436')

for s in [1, mp.mpf('0.5'), -mp.mpf('0.5'), 2]:
    target = mp.zeta(2, 1+mp.mpf(s)) if not hasattr(mp,'zeta') else None
    target = mp.zeta(2) if s==0 else None
    # Hurwitz zeta: zeta(2,1+s) = sum_{m>=1} 1/(m+s)^2 = psi'(1+s)
    target = mp.psi(1, 1+mp.mpf(s))
    vals = [(mp.mpf(1)/p, target - A_series_s(p, s)) for p in ps]
    coef = richardson(vals, 5)
    pred = [d1val, -(mp.mpf(s)+mp.mpf('0.5'))*mp.zeta(2), (6*s**2+6*s+1)*mp.zeta(3)/6, -mp.mpf('0.375')*s*(s+1)*(2*s+1)*mp.zeta(4)]
    print(f's={s}: extracted d1..d5:', [mp.nstr(c,10) for c in coef])
    print(f'      predicted d1..d4:', [mp.nstr(c,10) for c in pred])

# B variant
targetB = mp.mpf(3)/4*mp.zeta(2)
valsB = [(mp.mpf(1)/p, targetB - B_series(p)) for p in ps]
coefB = richardson(valsB, 5)
print('B: extracted d1..d5:', [mp.nstr(c,10) for c in coefB])
print('B: predicted d1..d3:', mp.nstr(d1val/2,10), ', 0,', mp.nstr(-13*mp.zeta(3)/24,10))
# PSLQ B's d4 against zeta4
d4B = coefB[3]
print('B d4 ~', mp.nstr(d4B,12), ' PSLQ vs zeta4:', mp.pslq([d4B, mp.zeta(4)], tol=mp.mpf('1e-8'), maxcoeff=10**6, maxsteps=10000))
'''
open('steps/step_16.py','w').write(code)
import subprocess
r = subprocess.run(['nice','-n','10','<project>/.venv/bin/python','steps/step_16.py'],capture_output=True,text=True,timeout=2400)
print(r.stdout); print(r.stderr[-800:] if r.returncode else '')
open('progress.log','a').write('step_16: variant expansions: A_s d2=-(s+1/2)z2, d3=(6s^2+6s+1)z3/6, d4=-3s(s+1)(2s+1)z4/8; B: d1/2, 0, -13z3/24\n')
