import mpmath as mp, subprocess, os
import cx6r_common as C
print("CX6R step_15: independent constant recomputation + small sec.7/8 checks")
mp.mp.dps = 100
a = mp.pi/(2*mp.sqrt(2))
x0 = 2*mp.exp(mp.pi/4)
K0 = 64/(mp.pi*a)*mp.exp(a/2+mp.pi/4)
T0 = 32/(3*mp.pi*a)*mp.exp(3*a/2+3*mp.pi/4)
Cstar = 2 + T0 + K0*(mp.exp(mp.pi/4)+mp.exp(x0)+mp.exp(mp.e*x0))
delta = 2-mp.sqrt(mp.pi)
conv = (3/(mp.e*delta))**3
total = mp.mpf('254.23') + mp.mpf('10600000')*conv
print("K0 =", mp.nstr(K0, 32), " (log: 70.0985853351116237515327912023)")
print("T0 =", mp.nstr(T0, 32), " (log: 170.657769858025927890213093913)")
print("C_* =", mp.nstr(Cstar, 32), " (log: 10578907.8235293971071051548801)")
print("C_* < 10,600,000:", Cstar < 10600000, " margin:", mp.nstr(10600000-Cstar, 10))
print("delta =", mp.nstr(delta, 32), "> 0:", delta > 0)
print("conversion =", mp.nstr(conv, 32), " (log: 114.096375458588480449897257230)")
print("total =", mp.nstr(total, 32), " (log: 1209421834.09103789276891092664)")
print("total < 1.21e9:", total < mp.mpf('1.21e9'), " margin:", mp.nstr(mp.mpf('1.21e9')-total, 10))
# H(1)=pi/4 and recurrence (sec.2 binomial identification)
print("H(1) =", mp.nstr(C.H(1), 20), " pi/4 =", mp.nstr(mp.pi/4, 20))
rec = max(abs(C.H(2*n+3)/C.H(2*n+1) - mp.mpf(2*n+1)/(2*n+4)) for n in range(10))
print("max recurrence deviation n<10:", mp.nstr(rec, 5))
# |H(r w)| <= sqrt(2)/r used in (20)
m = max(abs(C.H(r*C.w4))*r for r in [mp.mpf('0.5'),1,3,10,50,200])
print("max r|H(r w)| over test pts (bound sqrt(2)=", mp.nstr(mp.sqrt(2),10), "):", mp.nstr(m, 10))
# P_0 digamma identity (sec.8 / verify_rest.py)
def P0_int(s):
    return C.glint(lambda t: t**s/(1+t), [[0,0.5],[0.5,1]], 6)
def P0_dig(s):
    return (mp.digamma(1+s/2)-mp.digamma((s+1)/2))/2
d = max(abs(P0_int(mp.mpf('0.5')+r*C.w4)-P0_dig(mp.mpf('0.5')+r*C.w4)) for r in [mp.mpf('0.3'),2,9])
print("P_0 digamma identity max err:", mp.nstr(d, 5))
# D13s / D13u existence (supplied main-term certificates; contents not rerun)
base='<project>/research/08_theory/'
for d_ in ('D13s_mainterm_cert','D13u_mainterm_analytic','D13s','D13u'):
    p = os.path.join(base, d_)
    if os.path.isdir(p):
        print("dir exists:", d_, "->", sorted(os.listdir(p))[:8])
