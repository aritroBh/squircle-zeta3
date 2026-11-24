from mpmath import mp, mpf, mpc, log, pi, exp, fabs
mp.dps = 40
z0 = (-log(2)+pi*1j)/2
R = abs(z0)
def Q(t,s,u):
    return (1-u)**(2*t)*(s**(2*t)+(1-s)**(2*t)) + u**(2*t)
def dQ(t,s,u):
    return (2*(1-u)**(2*t)*(s**(2*t)*log(s)+(1-s)**(2*t)*log(1-s))
            + 2*log(1-u)*(1-u)**(2*t)*(s**(2*t)+(1-s)**(2*t)) + 2*log(u)*u**(2*t))
def zero(s,u,m=-1):
    x = -log(u)
    t = (-log(2)+(1j*pi*(2*m+1)))/(2*x)
    for _ in range(200):
        tn = t - Q(t,s,u)/dQ(t,s,u)
        if fabs(tn-t) < mpf(10)**(-32): return tn
        t = tn
    return t
# what step_06 sampled:
for u in [mpf('0.1'), mpf('0.01'), mpf('1e-3')]:
    x = -log(u); vals=[]
    for s in [mpf('0.5'), mpf('0.25'), mpf('0.1')]:
        t = zero(s,u); vals.append(float(fabs(1/t)*R/x))
    print(f"u={u}: |1/tau*|R/x at s=0.5,0.25,0.1:", [round(v,6) for v in vals])
# true min over s at those u (both branches, fine grid)
for u,x_ in [(mpf('0.1'),None),(mpf('0.01'),None),(mpf('0.05'),None)]:
    x = -log(u); mn=9; arg=None
    for i in range(2,99):
        s = mpf(i)/100
        for m in [-1,0]:
            t = zero(s,u,m); v = fabs(t)*x/R
            if v<mn: mn, arg = v, s
    print(f"u={u}: TRUE min_s |tau*|x/R = {mp.nstr(mn,8)} at s={arg}")
