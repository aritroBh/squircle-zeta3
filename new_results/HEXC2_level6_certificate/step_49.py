from flint import nmod_mat
M2 = nmod_mat([[1,2,3],[4,5,6],[7,8,13]], 1000000009)
rr, r = M2.rref(inplace=True)
print('rank', r)
print('rr[0,0]:', rr[0,0], type(rr[0,0]))
try:
    print('rr[0]:', rr[0])
except Exception as e:
    print('rr[0] error:', e)
print('entries:', [rr[0,j] for j in range(rr.ncols())])
print('tolist:', rr.tolist() if hasattr(rr,'tolist') else 'n/a')
