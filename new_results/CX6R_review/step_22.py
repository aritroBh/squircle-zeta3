import mpmath as mp
import cx6r_common as C
mp.mp.dps = 60
v1 = C.glint(lambda x: x**17, [[0,2]], 24)
v2 = C.glint(lambda x: mp.exp(-3*x), [[0,1],[1,5]], 40)
print("GL self-test x^17:", mp.nstr(abs(v1 - mp.mpf(2)**18/18), 5))
print("GL self-test exp:", mp.nstr(abs(v2 - (1-mp.exp(-15))/3), 5))
