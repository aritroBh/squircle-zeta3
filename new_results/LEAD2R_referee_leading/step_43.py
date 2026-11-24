import mpmath as mp
mp.mp.dps = 60
s = mp.mpf('0.99')
# known: int_0^1 x^{-s} (1+x) dx = 1/(1-s) + 1/(2-s)
v = mp.quad(lambda x: x**(-s)*(1+x), [0,1])
print("plain:", mp.nstr(v,12), "exact:", mp.nstr(1/(1-s)+1/(2-s),12), "err:", mp.nstr(abs(v-(1/(1-s)+1/(2-s))),3))
v2 = 1/(1-s)+mp.quad(lambda x: x**(-s)*((1+x)-1), [0,1])
print("sub1:", mp.nstr(v2,12), "err:", mp.nstr(abs(v2-(1/(1-s)+1/(2-s))),3))
# e^-r version
v3 = 1/(1-s) + mp.quad(lambda r: mp.exp(-(1-s)*r)*((1+mp.exp(-r))-1), [0,20,mp.inf])
print("e^-r sub1:", mp.nstr(v3,12), "err:", mp.nstr(abs(v3-(1/(1-s)+1/(2-s))),3))
v4 = mp.quad(lambda r: mp.exp(-(1-s)*r)*(1+mp.exp(-r)), [0, mp.inf])
print("e^-r plain:", mp.nstr(v4,12), "err:", mp.nstr(abs(v4-(1/(1-s)+1/(2-s))),3))
# F2 piece alone (no singularity): int_0^1 x^0.99 h dx computed two ways
h = lambda w: ((1+w)/mp.sqrt(1+w**2)-1)/w
F2a = mp.quad(lambda x: x**s*h(x), [0, mp.mpf('0.5'), 1])
F2b = mp.quad(lambda r: mp.exp(-(s+1)*r)*h(mp.exp(-r)), [0, mp.inf])
print("F2a:", mp.nstr(F2a,15), " F2b:", mp.nstr(F2b,15), " diff:", mp.nstr(abs(F2a-F2b),3))
