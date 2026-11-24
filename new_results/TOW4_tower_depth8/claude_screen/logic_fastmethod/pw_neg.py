import pw_adj as A
from fractions import Fraction as Fr
for (dl,s) in [(7,0),(7,2),(8,1),(8,3)]:
    r0=[A.residual(dl,s,Fr(p))[0].iszero() for p in range(4,10)]
    r1=[A.residual(dl,s,Fr(p),foldsign=-1)[0].iszero() for p in range(4,10)]
    r2=[A.residual(dl,s,Fr(p),rA0shift=1)[0].iszero() for p in range(4,10)]
    print((dl,s),'baseline all zero:',all(r0),'| fold sign flipped: zero at',sum(r1),'/6 | rA0 perturbed: zero at',sum(r2),'/6')
