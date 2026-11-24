
import mpmath as mp, time
mp.mp.dps = 50
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
h = mp.mpf('1e-12')
x, y = mp.mpf('0.5'), mp.mpf('0.7')
t=time.time()
fx = (li_convB(2,2,x+h,y) - li_convB(2,2,x-h,y))/(2*h)
print("fx", fx, time.time()-t)
# the second point y=1.4
x2, y2 = mp.mpf('0.31'), mp.mpf('1.4')
t=time.time()
v = li_convB(3,1,x2,y2)
print("Li31 at (0.31,1.4):", v, time.time()-t)
t=time.time()
v = li_convB(1,1,x2,y2)
print("Li11 at (0.31,1.4):", v, time.time()-t)
