import mpmath
from mpmath import mp, mpf, sqrt, nstr
mp.dps = 70
p = mpf(50); q = p/(p-1)
c1 = c2 = sqrt(2)/2
def F(z):
    return (-abs(c1-z)**(q-2)*(c1-z) - abs(c2-z)**(q-2)*(c2-z) + abs(z)**(q-2)*z)
for z in ['0','0.3','0.5','0.6','0.65','0.7','0.707','0.70710678118654752440','0.9','1.4']:
    print(z, nstr(F(mpf(z)),12))
