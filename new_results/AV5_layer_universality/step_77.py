import mpmath
from mpmath import mp, mpf, sqrt, pi, log, exp, cos, sin, tan, nstr
mp.dps = 70
def Fprof(th, lam):
    c, s = cos(th), sin(th)
    return c*log(lam*c) + (s/lam)*log(s) - (c+s/lam)*log(s+lam*c)
for lam in [mpf(2), mpf(3), mpf('0.5')]:
    d1 = 4*mpmath.quad(lambda th: Fprof(th, lam), [0, pi/4, pi/2])
    print(f"lam={lam}: d1 = 4 int F = {nstr(d1,45)}")
print("fit values: lam2 -3.0630013748005184342593242773606; lam3 -2.421476812773221095332596814398359; lam0.5 -6.1260027496010537996534024718262299")
# verify profile shape against measured values from ckpt
import json
rows = json.load(open('ckpt_lame_profile.json'))
lam = mpf(2)
print("\nprofile shape check (measured F(300) vs closed form):")
for th_s, v100, v300 in rows[::4]:
    th = mpf(th_s)
    print(f"  th={nstr(th,8)}: measured {nstr(mpf(v300),14)} vs closed {nstr(Fprof(th,lam),14)}")
