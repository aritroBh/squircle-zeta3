from mpmath import mp, mpf, pi, sin, cos, tan, log, quad, sqrt
import time
mp.dps = 60

one = mpf(1)
def H1u(u):
    s = sin(u); c = cos(u)
    if abs(c) < mpf('1e-30'):  # u near pi/2: limit forms
        # H1(sin u)*cos u:  (2 s^2 ln s - (s^2-1) ln(1-s^2))/(4(1-s^2)^(3/2)) * c
        # (1-s^2)=c^2: = (2 s^2 ln s - (s^2-1)*2 ln c)/(4 c^3) * c = (2s^2 ln s + 2 c^2 ln c... careful
        pass
    s2 = s*s
    val = (2*s2*log(s) - (s2-1)*log(1-s2))/(4*(1-s2)**mpf('1.5'))
    return val*c

def H2u(u):
    s = sin(u); c = cos(u); s2 = s*s; c2 = 1-s2
    ls = log(s); lc2 = log(c2)
    A = 2*s2*ls - (s2-1)*lc2   # = 2 s2 ls + c2 lc2
    num = 8*s2*s2*ls**2 - 8*s2*(s2-1)*(ls-1)*ls - 4*c2**2*lc2 + A**2
    return num/(16*(c2)**mpf('2.5'))*c

def H3u(u):
    s = sin(u); c = cos(u); s2 = s*s; c2 = 1-s2
    ls = log(s); lc2 = log(c2)
    A = 2*s2*ls - (s2-1)*lc2
    B = -2*s2*s2*ls**2 + 2*s2*(s2-1)*(ls-1)*ls + c2**2*lc2
    num = (64*s2**3*ls**3 - 48*s2*s2*(s2-1)*(2*ls-1)*ls**2 - 48*s2*(s2-1)**2*(ls-1)*ls
           + 32*s2*(s2-1)**2*ls**3 - 24*c2**3*lc2 + A**3 - 12*A*B)
    return num/(64*(c2)**mpf('3.5'))*c

# The above have 0/0 at u=pi/2 (c->0); check behavior. Limits:
# H1u ~ (2*1*ln(1) + c2*ln c2)/(4 c2^1.5) *c = (2 c2 ln c)/(4 c^3)*c = ln(c)/2 -> -inf? log singularity, integrable
# Let's just test numerically near pi/2
for u in [1.5, 1.57, 1.5707]:
    print(u, mp.nstr(H1u(u),10), mp.nstr(H2u(u),10), mp.nstr(H3u(u),10))
