"""Step 2 (exact): Mobius maps phi with phi({0,1,-1,inf,b^2}) inside S={0,+-1,+-a,+-b,inf}.
phi is fixed by the images of 0,1,inf (336 ordered triples)."""
import itertools
from qs2 import Q, A, B, R2
INF=None
S=[Q(0),Q(1),Q(-1),A,-A,B,-B,INF]
nm={0:'0',1:'1',2:'-1',3:'a',4:'-a',5:'b',6:'-b',7:'inf'}
def mob(p,q,r):
    """phi with phi(0)=p, phi(1)=q, phi(inf)=r; returns function on Q or None(=inf)."""
    # phi = M^-1 o psi where psi(z)=CR sending p,q,r -> 0,1,inf: psi(z)=(z-p)(q-r)/((z-r)(q-p))
    # solve psi(w)=t for w.  Write as matrix.
    import fractions
    def mat_psi():
        # psi(z) = (z-p)(q-r) / ((z-r)(q-p)) ; handle inf among p,q,r
        if r is INF: return [[Q(1),-p],[Q(0),q-p]]
        if p is INF: return [[Q(0),q-r],[Q(1),-r]]
        if q is INF: return [[Q(1),-p],[Q(1),-r]]
        return [[q-r,-(q-r)*p],[q-p,-(q-p)*r]]
    M=mat_psi()
    (a,b),(c,d)=M
    # inverse: w = (d t - b)/(-c t + a)
    def f(t):
        if t is INF:
            num,den=d,-c
        else:
            num,den=d*t-b,-c*t+a
        if den.iszero(): return INF
        return num/den
    return f
def key(x): return 'inf' if x is INF else repr(x)
found=[]
for i,j,k in itertools.permutations(range(8),3):
    f=mob(S[i],S[j],S[k])
    assert key(f(Q(0)))==key(S[i]) and key(f(Q(1)))==key(S[j]) and key(f(INF))==key(S[k])
    m1=f(Q(-1)); z=f(B*B)
    keys=[key(s) for s in S]
    if key(m1) in keys and key(z) in keys:
        found.append((nm[i],nm[j],nm[k],nm[keys.index(key(m1))],nm[keys.index(key(z))]))
print('phi(0),phi(1),phi(inf),phi(-1),phi(b^2):')
for t in found: print(t)
print('count',len(found))
