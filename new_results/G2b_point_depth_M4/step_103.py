
import mpmath as mp
mp.mp.dps = 120

def reli2(z):
    z = mp.mpf(z)
    if z == 0: return mp.mpf(0)
    if abs(z) <= mp.mpf('0.55'):
        tot = mp.mpf(0); zk = z; k = 1
        while True:
            t = zk / (k * k)
            tot += t
            if abs(t) < mp.mpf(10) ** (-mp.mp.dps + 10): break
            zk *= z; k += 1
        return tot
    if z > 1:
        return -reli2(1 / z) + mp.pi**2 / 3 - mp.log(z)**2 / 2
    if z > 0:  # (0.55, 1]
        return mp.pi**2 / 6 - mp.log(z) * mp.log(1 - z) - reli2(1 - z)
    # z < 0: Landen
    w = z / (z - 1)
    return -mp.log(1 - z)**2 / 2 - reli2(w)

for z in [mp.mpf('0.2'), mp.mpf('-0.2'), mp.mpf('0.8'), mp.mpf('1.5'), mp.mpf('-2'), mp.mpf('5'),
          mp.mpf('-0.9'), mp.mpf('0.999'), mp.mpf('-0.5'), mp.mpf('-100'), mp.mpf('0.55')]:
    a = reli2(z); b = mp.polylog(2, z).real
    print(z, mp.nstr(abs(a-b), 3))
