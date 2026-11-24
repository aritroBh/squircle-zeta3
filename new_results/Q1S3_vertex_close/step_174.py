import mpmath as mp
mp.mp.dps=100
x=mp.mpf('0.22799665392074099901')
print(mp.pslq([x, x], tol=mp.mpf('1e-80'), maxcoeff=100))
print(mp.pslq([x, 2*x], tol=mp.mpf('1e-80'), maxcoeff=100))
# with a bigger vector
import mpmath
vec=[x, x, mp.pi, mp.e]
print(mp.pslq(vec, tol=mp.mpf('1e-80'), maxcoeff=100))
