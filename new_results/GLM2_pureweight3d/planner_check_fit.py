import mpmath as mp
mp.mp.dps=50
pts=[]
for line in open('sp_values.txt'):
    if line.startswith('#') or not line.strip(): continue
    p,S=line.split(); pts.append((mp.mpf(p),mp.mpf(S)))
r3=mp.sqrt(3); pi=mp.pi; z3=mp.zeta(3); z5=mp.zeta(5)
claim=lambda t: -(pi**2/9)*t**2-(16*z3/9)*t**3-(37*pi**4/810)*t**4-(160*z5/27)*t**5
print('S(1)-4sqrt3 =', mp.nstr([S for p,S in pts if p==1][0]-4*r3,3))
rows=[];rhs=[]
for p,S in pts:
    if p==1: continue
    t=1-1/p; F=mp.log((2-p)*S/(4*r3)); R=F-claim(t)
    print(f'p={mp.nstr(p,4):6s} tau={mp.nstr(t,6):10s} residual={mp.nstr(R,6):12s} residual/tau^6={mp.nstr(R/t**6,6)}')
# free fit: F = sum_{n=1..8} c_n tau^n
for p,S in pts:
    if p==1: continue
    t=1-1/p; rows.append([t**n for n in range(1,9)]); rhs.append(mp.log((2-p)*S/(4*r3)))
c=mp.lu_solve(mp.matrix(rows),mp.matrix(rhs))
for n,ref in [(1,0),(2,-pi**2/9),(3,-16*z3/9),(4,-37*pi**4/810),(5,-160*z5/27)]:
    print(f'c{n} = {mp.nstr(c[n-1],12):18s} claim {mp.nstr(ref,12):18s} diff {mp.nstr(c[n-1]-ref,3)}')
