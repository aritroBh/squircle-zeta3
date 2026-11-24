import mpmath, json
from mpmath import mp, mpf, pi, nstr
mp.dps = 240
C = json.load(open('ckpt_constants.json'))
def gv(k): return mpf(C[k])
lam, L, c = gv('lam'), gv('L'), gv('c8')
l1mc = mp.log(1-c)
res = json.load(open('ckpt_Mk_cornerintegral.json'))
M0 = mpf(res['3pi/4'][0])
test = 2*M0 - 3*lam + 2*L - 3*c*lam + 2*c*L - 4*c*l1mc
print("2*M0 -3lam +2L -3c*lam +2c*L -4c*ln(1-c) =", nstr(test, 10))
# so the relation exists: 2*M0 = 3*lam - 2*L + 3*c*lam - 2*c*L + 4*c*ln(1-c)
