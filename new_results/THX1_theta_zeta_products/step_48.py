
import flint
A = flint.nmod_mat(2, 3, 7)
A[0,0]=1; A[0,1]=2; A[1,0]=3; A[1,2]=1
r = A.rref()
print(type(r), r if not hasattr(r,'__len__') else len(r))
if isinstance(r, tuple):
    for x in r: print(type(x))
