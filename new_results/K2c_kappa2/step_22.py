
import mpmath
from mpmath import mp, mpf, mpc, sqrt, log, exp, pi, gamma

def run(dps, a, N, m):
    mp.dps = dps
    u4 = exp(mpc(0,1)*mp.pi/4)
    Z = 1j*a; aa = abs(Z); s2 = mp.sqrt(2)
    f0 = lambda t: 1/s2 - (N/2)*(s2*aa*t + 2*aa*aa)/(t**3 + s2*aa*t*t + aa*aa*t)
    ts = mp.findroot(f0, mp.sqrt(N*aa))
    def ft(t):
        x = t*u4
        return exp(-x)*(1+Z/x)**(-N)*(x+Z)**(m+mpf(3)/2)/x * u4
    w = ts/6
    pieces = [[0, ts*mpf('0.3')],[ts*mpf('0.3'), ts*mpf('0.6')],[ts*mpf('0.6'), ts-w],[ts-w, ts-w/3],[ts-w/3, ts],[ts, ts+w/3],[ts+w/3, ts+w],[ts+w, ts*mpf('1.6')],[ts*mpf('1.6'), ts*mpf('3')],[ts*mpf('3'), ts*mpf('10')],[ts*mpf('10'), mp.inf]]
    tot = 0
    vals = []
    for p in pieces:
        v = mp.quadgl(ft, p)
        vals.append(v); tot += v
    return ts, vals, tot

a = 2*mp.pi/3; N = 5600
for m in [0,1,2]:
    ts, v70, t70 = run(70, a, N, m)
    ts, v110, t110 = run(120, a, N, m)
    print('m=%d ts=%s' % (m, mp.nstr(ts, 10)))
    for i,(x,y) in enumerate(zip(v70, v110)):
        d = abs(x-y)
        if d > abs(t110)*mpf('1e-25'):
            print('  piece %d: dps70=%s dps120=%s' % (i, mp.nstr(x,8), mp.nstr(y,8)))
    print('  total70  =', mp.nstr(t70, 15))
    print('  total120 =', mp.nstr(t110, 15))
    print('  diff =', mp.nstr(abs(t70-t110), 4))
