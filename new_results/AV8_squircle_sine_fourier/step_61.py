from fractions import Fraction as Fr
# Values are q0 + q1*pi + q2*pi*ln2 ; represent as tuple (q0,q1,q2) of Fractions
def A(m):  # = pi/2
    return (Fr(0), Fr(1,2), Fr(0))
def B(m):  # m even >=2: -pi/2 ln2 + (-1)^{m/2} pi/(2m)
    k = m//2
    return (Fr(0), (-1)**k*Fr(1,2*m), Fr(-1,2))
def K(r):  # r odd
    if r == 1: return (Fr(0), Fr(1,8), Fr(-1,4))
    return (Fr(0), -Fr(1,2*(r*r-1)), Fr(0))
def Kt(m):  # m=2k: 2*sum_{j=1}^k (-1)^{j-1} K(m-2j+1)
    k = m//2
    acc = [Fr(0),Fr(0),Fr(0)]
    for j in range(1,k+1):
        sgn = (-1)**(j-1)
        q = K(m-2*j+1)
        for i in range(3): acc[i] += 2*sgn*q[i]
    return tuple(acc)
def Tval(m):  # m=2k: -2k Kt - pi/2 + pi/2 ln2 - pi/(4k)
    k = m//2
    kt = Kt(m)
    return (-2*k*kt[0], -2*k*kt[1] - Fr(1,2) - Fr(1,4*k), -2*k*kt[2] + Fr(1,2))
def add(*vs):
    return tuple(sum(v[i] for v in vs) for i in range(3))
def bn1(n):
    # (1/pi) * sum_{m in {n-1,n+1}} (A+B+T)(m)/m ; n odd >=3
    tot = [Fr(0),Fr(0),Fr(0)]
    for m in [n-1,n+1]:
        v = add(A(m),B(m),Tval(m))
        for i in range(3): tot[i] += v[i]/m
    # divide by pi: q0+q1*pi+q2*pi*ln2 -> q1 + q2*ln2
    return tot[1], tot[2]   # b_n' = q1 + q2*ln2

for n in [3,5,7,9,11,13,15,17,19,21]:
    q1,q2 = bn1(n)
    print(n, "b_n' =", q1, "+", q2, "*ln2", "  =", float(q1+q2*Fr(6931471805599453,10**16)))
