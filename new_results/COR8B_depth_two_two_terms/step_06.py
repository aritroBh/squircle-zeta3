
import mpmath as mp, time
mp.mp.dps = 50
t = time.time()
print(mp.polylog(2, mp.mpf('0.5')), time.time()-t)
t = time.time()
# test li_convB quickly
def li_convB(a, b, x, y):
    x = mp.mpf(x); y = mp.mpf(y)
    total = mp.mpf(0); m = 0; partial = mp.mpf(0); xp = mp.mpf(1)
    Li_a_x = mp.polylog(a, x)
    tol = mp.mpf(10)**(-mp.mp.dps+8)
    while True:
        m += 1; xp *= x
        partial += xp / mp.power(m, a)
        term = mp.power(y, m)/mp.power(m, b) * (Li_a_x - partial)
        total += term
        if m > 40 and abs(term) < tol: break
    return total
print(li_convB(2,2,mp.mpf('0.5'),mp.mpf('0.7')), time.time()-t)
