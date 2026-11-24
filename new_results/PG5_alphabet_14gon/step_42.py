import mpmath as mp
mp.mp.dps = 50
m = 7
nn = [mp.matrix([mp.cos(i*mp.pi/m), mp.sin(i*mp.pi/m)]) for i in range(m)]
# K1: sum_i |<n_i,x>| <= 1. Edge in sector (-pi/14, pi/14): sum = <u,x>, u = n0+n1+n2+n3-n4-n5-n6
u = nn[0]+nn[1]+nn[2]+nn[3]-nn[4]-nn[5]-nn[6]
print("u =", u, " |u| =", mp.norm(u))
# apothem h = 1/|u|, edge x1 = 1/U with U=u[0]
h = 1/mp.norm(u)
print("h =", mp.nstr(h,30), " tan(pi/14)=L?", mp.nstr(mp.tan(mp.pi/14),30))
V0 = mp.matrix([mp.mpf('0.5'), -mp.mpf('0.5')*mp.tan(mp.pi/14)])
V1 = mp.matrix([mp.mpf('0.5'),  mp.mpf('0.5')*mp.tan(mp.pi/14)])
L = mp.norm(V1-V0)
print("L =", mp.nstr(L,30), " h^2/L^2 =", mp.nstr(h*h/(L*L),30))
# term values a_i(t) = |<n_i, x(t)>| with positive sign
def x(t): return V0 + t*(V1-V0)
# signs in the sector: i=0,1,2,3 positive; 4,5,6 negative
sgn = {0:1,1:1,2:1,3:1,4:-1,5:-1,6:-1}
def a(i,t): return sgn[i]*(nn[i][0]*x(t)[0] + nn[i][1]*x(t)[1])
# check sum = 1 and linearity
for tt in ['0','0.3','0.7','1']:
    tt = mp.mpf(tt)
    print("sum a_i at", tt, "=", mp.nstr(sum(a(i,tt) for i in range(7)),20))
# term data: a_i(t) = B_i(1-x_i t) for nonlocals; locals i=3 (A0 t), i=4 (A1(1-t))
A0 = a(3,1)  # a3(t)=A0*t since a3(0)=0
A1 = a(4,0)
print("A0 =", mp.nstr(A0,30)); print("A1 =", mp.nstr(A1,30))
for i in [1,2,5,6]:
    B = a(i,0); slope = a(i,1)-a(i,0); xi = -slope/B
    print(f"i={i}: B={mp.nstr(B,25)} x={mp.nstr(xi,25)} 1-x={mp.nstr(1-xi,25)} p=1/x={mp.nstr(1/xi,25)}")
print("constant term i=0:", mp.nstr(a(0,0),20))
# Landen pairs: x6 vs x1/(x1-1)
x1 = -(a(1,1)-a(1,0))/a(1,0); x6 = -(a(6,1)-a(6,0))/a(6,0)
x2 = -(a(2,1)-a(2,0))/a(2,0); x5 = -(a(5,1)-a(5,0))/a(5,0)
print("x6 - x1/(x1-1):", mp.nstr(x6 - x1/(x1-1),5))
print("x5 - x2/(x2-1):", mp.nstr(x5 - x2/(x2-1),5))
# check P'(1)... apothem vs m=6 style: h should be 1/2
