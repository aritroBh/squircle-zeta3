from flint import nmod_mat
M = nmod_mat([[1,2,3],[4,5,6],[7,8,10]], 1000000009)
out = M.rref(inplace=True)
print(type(out), out)
M2 = nmod_mat([[1,2,3],[4,5,6],[7,8,10]], 1000000009)
rr, r = M2.rref(inplace=True)
print(type(rr), rr if not hasattr(rr,'__len__') else 'len', r)
print(repr(rr))
