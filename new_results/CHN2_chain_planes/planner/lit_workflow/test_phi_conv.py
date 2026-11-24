import sys; sys.dont_write_bytecode=True
sys.path.insert(0,'.')
import random
from fractions import Fraction as F
import phi
random.seed(1)
n=6
def rnd(): return [[F(random.randint(-3,3)) for _ in range(n)] for _ in range(n)]
def mm(A,B): return [[sum(A[i][k]*B[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
def add(A,B,c=1): return [[A[i][j]+c*B[i][j] for j in range(n)] for i in range(n)]
def com(A,B): return add(mm(A,B),mm(B,A),-1)
e0=rnd(); e1=rnd()
out=phi.phi_first_column(e0,e1,5)
# weight 2: zeta2 [e0,e1]; weight 3: zeta3 [e0-e1,[e0,e1]]; weight 4 per BSST (11)
c01=com(e0,e1)
W2=c01
W3=com(add(e0,e1,-1),c01)
W4=add(add(add(com(e0,com(e0,com(e0,e1))), com(e1,com(e0,com(e1,e0))),F(1,4)), com(e1,com(e1,com(e1,e0))),-1), mm(c01,c01),F(5,4))
for W,Mx in ((2,W2),(3,W3),(4,W4)):
    col=[Mx[i][0] for i in range(n)]
    got=[out[W][i][0] for i in range(n)]
    print(W, 'basis=',phi.E.table(W).monomials,'match:', col==got)
