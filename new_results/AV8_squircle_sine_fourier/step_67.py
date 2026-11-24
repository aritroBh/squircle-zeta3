from fractions import Fraction as Fr
# verified Fraction code pieces
def A(m): return (Fr(0), Fr(1,2), Fr(0))
def B(m):
    k = m//2
    return (Fr(0), (-1)**k*Fr(1,2*m), Fr(-1,2))
def K(r):
    if r == 1: return (Fr(0), Fr(1,8), Fr(-1,4))
    return (Fr(0), -Fr(1,2*(r*r-1)), Fr(0))
def Kt(m):
    k = m//2
    acc = [Fr(0),Fr(0),Fr(0)]
    for j in range(1,k+1):
        sgn = (-1)**(j-1)
        q = K(m-2*j+1)
        for i in range(3): acc[i] += 2*sgn*q[i]
    return tuple(acc)
def Tval(m):
    k = m//2
    kt = Kt(m)
    return (-2*k*kt[0], -2*k*kt[1] - Fr(1,2) - Fr(1,4*k), -2*k*kt[2] + Fr(1,2))
def add(*vs): return tuple(sum(v[i] for v in vs) for i in range(3))
for mm in [4,6]:
    v = add(A(mm),B(mm),Tval(mm))
    print(mm, "A+B+T (q_pi, q_piL):", v[1], v[2], " /m:", v[1]/mm, v[2]/mm)
