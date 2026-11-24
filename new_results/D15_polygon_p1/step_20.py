import mpmath as mp
mp.mp.dps = 50
m = 2; A = mp.pi/m
def S(th, d):
    p = 1+d
    return mp.fsum([mp.sin(th + i*A)**p for i in range(m)])
def logrho(th, d):
    p = 1+d
    return -mp.log(S(th,d))/p
def Q_num(th, d):
    h2 = mp.mpf('1e-12')
    return (logrho(th-2*h2,d) - 8*logrho(th-h2,d) + 8*logrho(th+h2,d) - logrho(th+2*h2,d))/(12*h2)
th = mp.mpf('0.3')
h = mp.mpf('1e-9')
q = lambda dd: Q_num(th, dd)
q0 = q(0)
q1 = (q(-2*h) - 8*q(-h) + 8*q(h) - q(2*h))/(12*h)
print("Q0 num:", mp.nstr(q0,25))
print("Q1 num:", mp.nstr(q1,25))

# formula pieces
sv = [mp.sin(th+i*A) for i in range(m)]; cv = [mp.cos(th+i*A) for i in range(m)]
ls = [mp.log(v) for v in sv]
T = [mp.fsum([sv[i]*ls[i]**k for i in range(m)]) for k in range(3)]
Tp = [mp.fsum(cv)] + [mp.fsum([cv[i]*(ls[i]**k + k*ls[i]**(k-1)) for i in range(m)]) for k in range(1,3)]
S0 = T[0]; S0p = Tp[0]
U1 = T[1]/S0
U1p = Tp[1]/S0 - T[1]*S0p/S0**2
Q0f = S0p/S0
Q1f = mp.log(S0)*0  # placeholder
Q1f = (S0p/S0) - U1p   # R1' = (log S0)' - U1' = Q0 - U1p
print("Q0 formula:", mp.nstr(Q0f,25))
print("Q1 formula:", mp.nstr(Q1f,25))
print("U1p formula:", mp.nstr(U1p,25))
# numeric U1p:
def U1_num(th):
    svv = [mp.sin(th+i*A) for i in range(m)]
    lsv = [mp.log(v) for v in svv]
    return mp.fsum([svv[i]*lsv[i] for i in range(m)])/mp.fsum(svv)
h2 = mp.mpf('1e-12')
U1p_num = (U1_num(th-2*h2) - 8*U1_num(th-h2) + 8*U1_num(th+h2) - U1_num(th+2*h2))/(12*h2)
print("U1p numeric:", mp.nstr(U1p_num,25))
