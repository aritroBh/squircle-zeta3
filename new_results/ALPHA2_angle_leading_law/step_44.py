from mpmath import mp, mpf, quad, sqrt, nstr, pi, sin, cos, gamma, legenp
mp.dps = 50
def Qf(w, ca): return sqrt(1 - 2*w*ca + w**2)
def Eint(s, ca): return quad(lambda w: (w**(1+s) + w**(-s))/Qf(w,ca)**3, [0, 1])
def Ecand2(s, ca, alpha):
    return -pi*legenp(-1-s, 1, -ca, type=2)/(sin(alpha)*sin(pi*s))
for alpha in [pi/2, pi/3, 2*pi/3, pi/6]:
    ca = cos(alpha)
    for s in [mpf('0.3'), mpf('-0.2'), mpf('0.55')]:
        a = Eint(s, ca); b = Ecand2(s, ca, alpha)
        print(f"alpha={nstr(alpha,6)}, s={nstr(s,5)}: |E - [-pi P1/(sa sin ps)]|={nstr(abs(a-b),3)}")
import subprocess
print(subprocess.run(['tail','-3','<project>/research/08_theory/ALPHA2_angle_leading_law/kern_old_ext.out'], capture_output=True, text=True).stdout)
