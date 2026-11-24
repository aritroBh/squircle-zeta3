from flint import fmpz_mat
A = fmpz_mat([[2,0],[0,3]])
b = fmpz_mat([[4],[9]])
out = fmpz_mat.solve(A, b)
print(type(out), out)
help(fmpz_mat.solve)
