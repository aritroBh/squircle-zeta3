from mpmath import mp, mpf, quad, gamma, hyp3f2, nstr
mp.dps = 30
A = [3, 1, 2, 5]; u = mpf(1)/7
u1, u2, u3, u4 = [a*u for a in A]; p1, p2, p3, p4 = [1/x for x in (u1, u2, u3, u4)]
# x-coordinates: x2 in [0,1], x3 in [0, (1-x2^p2)^(1/p3)]
d = quad(lambda a: (1-a**p2)**u1*quad(lambda b: (1-b**p3)**u4, [0, (1-a**p2)**(1/p3)]), [0, mpf(1)/2, mpf('0.9'), 1])
F = gamma(1+u2)*gamma(1+u1+u3)/gamma(1+u1+u2+u3)*hyp3f2(-u4, u3, 1+u1+u3, 1+u3, 1+u1+u2+u3, 1)
print(nstr(d, 25), nstr(F, 25), nstr(d-F, 3))
