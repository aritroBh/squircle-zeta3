import mpmath as mp
mp.mp.dps = 60
def L_direct(chimap, N=4000):
    return mp.nsum(lambda n: chimap(n%8)/n**5 if n%8 else 0, [1, N])
# convention A: chi(1)=chi(3)=1, chi(5)=chi(7)=-1
A = L_direct(lambda r: 1 if r in (1,3) else -1)
# convention B: chi(1)=1, chi(3)=-1, chi(5)=-1, chi(7)=1
B = L_direct(lambda r: 1 if r in (1,7) else -1)
print('A (1,3)+:', mp.nstr(A, 25))
print('B (1,7)+:', mp.nstr(B, 25))
print('planner : 0.9956338563129674563422177171275')
# also via hurwitz for B
LB = (mp.zeta(5, mp.mpf(1)/8) - mp.zeta(5, mp.mpf(3)/8) - mp.zeta(5, mp.mpf(5)/8) + mp.zeta(5, mp.mpf(7)/8))/8**5
print('B hurwitz:', mp.nstr(LB, 25))
