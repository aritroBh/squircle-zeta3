import subprocess
r = subprocess.run(['<project>/.venv/bin/python','-c',
'''
import flint
# tall matrix 5x3, rank 3
M = flint.fmpz_mat([[1,2,3],[4,5,6],[7,8,10],[1,0,1],[2,2,2]])
H = M.hnf()
for i in range(H.nrows()):
    print([int(H[i,j]) for j in range(H.ncols())])
print("rank:", M.rank())
'''], capture_output=True, text=True, cwd='.')
print(r.stdout, r.stderr)
