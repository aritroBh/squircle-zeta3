import sys, time
sys.path.insert(0, '.')
import mpmath as mp
import cor7_lib as L

mp.mp.dps = 1150
c = mp.sqrt(2) - 1
t0 = time.time()
M3, W4, W3, W2 = L.M3_weights(c)
print(f"HPL at dps 1150: {time.time()-t0:.0f}s")
for nm, v in [("M3", M3), ("W4", W4), ("W3", W3), ("W2", W2)]:
    open(f"b21_{nm}_hpl_1150.txt", "w").write(mp.nstr(v, 1145))
    print(nm, "=", mp.nstr(v, 40))
alpha = 2*mp.acos(c)
print("alpha =", mp.nstr(alpha, 40))
print("alpha/pi =", mp.nstr(alpha/mp.pi, 30))
print("alpha deg =", mp.nstr(alpha*180/mp.pi, 20))
