# endpoint behaviour of N_{(i,j),a} at v=a: log t vanishes linearly in (a-v); log(1-t) has a log singularity
import mpmath as mp
mp.mp.dps = 60
s2 = mp.sqrt(2); a = 1 + s2; b = s2 - 1
for d in [mp.mpf(10) ** -4, mp.mpf(10) ** -8, mp.mpf(10) ** -16, mp.mpf(10) ** -30]:
    v = a - d
    omt = d * (v + b) / (4 * v)
    logt = mp.log1p(-omt)
    print('a-v=%s  log t/(v-a) = %s   log(1-t)=%s' % (mp.nstr(d, 3), mp.nstr(logt / (v - a), 15), mp.nstr(mp.log(omt), 8)))
print('predicted limit: -(a+b)/(4a)*(-1) = (a+b)/(4a) =', mp.nstr((a + b) / (4 * a), 15), ' (log t ~ -(a-v)(a+b)/(4a))')
# integrand magnitude of log^i t log^j(1-t)/(v-a) near the endpoint, for the three keys
for (i, j) in [(5, 0), (4, 1), (3, 2)]:
    for d in [mp.mpf(10) ** -6, mp.mpf(10) ** -12]:
        v = a - d
        omt = d * (v + b) / (4 * v)
        val = mp.log1p(-omt) ** i * mp.log(omt) ** j / (v - a)
        print('key', (i, j), 'a-v=%s' % mp.nstr(d, 2), 'integrand =', mp.nstr(val, 6))
