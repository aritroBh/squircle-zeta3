import time, sys
sys.path.insert(0, '.')
import mpmath as mp
import cor7_lib as L

mp.mp.dps = 1150
c3 = mp.mpf(1)/3
t0 = time.time()
M3, W4, W3, W2 = L.M3_weights(c3)
print(f"HPL at dps 1150: {time.time()-t0:.0f}s")
for nm, v in [("M3", M3), ("W4", W4), ("W3", W3), ("W2", W2)]:
    open(f"c13_{nm}_hpl_1150.txt", "w").write(mp.nstr(v, 1145))
    print(nm, "=", mp.nstr(v, 40))
print("W4+W3+W2-M3 =", mp.nstr(abs(W4+W3+W2-M3), 3))
