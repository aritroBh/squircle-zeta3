import mpmath as mp
mp.mp.dps = 50
# test 1: plain singular integral int_0^1 x^-0.99 dx = 100
v1 = mp.quad(lambda x: x**(-mp.mpf('0.99')), [0,1])
print("int x^-0.99 =", mp.nstr(v1,10), "err", mp.nstr(abs(v1-100),3))
v1b = mp.quad(lambda x: x**(-mp.mpf('0.99')), [0, mp.mpf('1e-6'), 1])
print("with split:", mp.nstr(v1b,10), "err", mp.nstr(abs(v1b-100),3))
# WW at s=0.99 with split points
q = lambda w: (1+w)/(1+w**2)**mp.mpf('1.5')
s = mp.mpf('0.99')
W1 = mp.quad(lambda x: (x**s+x**(-s))*q(x), [0, mp.mpf('1e-8'), mp.mpf('0.5'), 1])
rW = mp.gamma(1+s/2)*mp.gamma(mp.mpf('0.5')-s/2)/mp.sqrt(mp.pi) + mp.gamma(1-s/2)*mp.gamma(mp.mpf('0.5')+s/2)/mp.sqrt(mp.pi)
print("WW with splits:", mp.nstr(W1, 12), "rhs:", mp.nstr(rW,12), "err:", mp.nstr(abs(W1-rW),3))
W2 = mp.quad(lambda x: x**(-s)*q(x), [0, mp.mpf('1e-8'), 1]) + mp.quad(lambda x: x**s*q(x), [0,1])
print("WW two-piece:", mp.nstr(W2, 12), "err:", mp.nstr(abs(W2-rW),3))
