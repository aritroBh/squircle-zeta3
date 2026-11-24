import mpmath as mp
mp.mp.dps = 60

def Li11(a,b):  # sum_{n>m>=1} a^n b^m/(n m)
    return mp.nsum(lambda n: a**n * mp.nsum(lambda m: b**m/m, [1, n-1], tol=mp.mpf('1e-70'))/n,
                   [2, mp.inf], tol=mp.mpf('1e-55'))

# octagon kappa values
sq2 = mp.sqrt(2)
kap1 = sq2-1; kap2 = 1/sq2-1
x1 = -kap1; x2 = -kap2
print("kappa1 =", mp.nstr(kap1,30), " kappa2 =", mp.nstr(kap2,30))
print("x1 = 1-sqrt2 =", mp.nstr(x1,30), " x2 = 1-1/sqrt2 =", mp.nstr(x2,30))
print("Landen check x2 vs x1/(x1-1):", mp.nstr(x2 - x1/(x1-1),3))

# direct integral I111
I111 = mp.quad(lambda t: mp.log(1+kap1*t)*mp.log(1+kap2*t), [0,1])
print("\nI111 direct =", mp.nstr(I111,40))

# my reduction: I111 = 2*Li11(x2,x1/x2) - 2*U
# U = (-log(1-x1))*(-log(1-x2)-x2)/x2 - (1/x2)*(Li11(x1,x2/x1) - x2*(-log(1-x1)-x1)/x1)
l1 = -mp.log(1-x1); l2_ = -mp.log(1-x2)
A = Li11(x2, x1/x2); Bv = Li11(x1, x2/x1)
U = l1*(l2_-x2)/x2 - (1/x2)*(Bv - x2*(l1-x1)/x1)
I111_formula = 2*A - 2*U
print("I111 formula =", mp.nstr(I111_formula,40))
print("diff:", mp.nstr(I111-I111_formula,3))
print("\nLi11(x2,x1/x2) =", mp.nstr(A,40))
print("Li11(x1,x2/x1) =", mp.nstr(Bv,40))
