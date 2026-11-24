# pointwise check of the claimed equivalence  closing(delta)  <=>  a_delta = feed(U), U(dp)=sum_s (-1)^s F(dp,s)
import random
import pw
from pw import G, Fu, TAU1, TAU2, I_, rand_pt, a_val
rng=random.Random(99)
def U(dp,D,M):
    tot=G(0)
    for s in range(0, dp//2+1): tot = tot + ((-1)**s)*Fu(dp,s,D,M)
    return tot
def feedU(delta,D,M):
    tot=G(0)
    for dp in range(0,delta+1):
        tot = tot + I_*(2*D-1)/(2*D-dp-1)*(TAU1[delta-dp]*U(dp,D+1,M+1)-TAU2[delta-dp]*U(dp,D+1,M))
    return tot
for delta in range(0,9):
    ok=True; clos=True
    for _ in range(5):
        D,M=rand_pt(rng)
        ok &= (a_val(delta,D,M)-feedU(delta,D,M)).iszero()
        clos &= Fu(delta, delta//2+1, D, M).iszero()
    print(f'delta={delta}: closing F(delta,top+1)=0 : {clos};   a_delta == feed(U) : {ok}')
# and a negative control: break a_7 by adding a d-dependent perturbation -> both closing and master should fail
