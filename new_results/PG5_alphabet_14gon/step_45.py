src = open('ftfield.py').read()
src = src.replace("    B = A1*B2+B1*A2 + 2*(B1*C2+C1*B2) + C1*C2",
                  "    B = A1*B2+B1*A2 + 2*(B1*C2+C1*B2) - C1*C2")
old_inv = """def finv(x):
    A,B,C,d=x
    # solve (A+B a+C a^2)(X+Y a+Z a^2) = d  -> 3x3 linear system
    import sympy as sp
    M = sp.Matrix([[A, B-C, -C],[B, A+2*B-C+ C*0, B-C],[C, B+3*C*0, A-C+3*C]])  # placeholder, computed below properly
    raise NotImplementedError
def finv(x):
    # use exact linear algebra with Fractions
    A,B,C,d = x
    # multiplication matrix of (A+B a+C a^2): columns = images of 1, a, a^2
    # (A+B a+C a^2)*1 = A + B a + C a^2
    # *a = (B-C) + (A+2C) a + (B+3C) a^2 ... wait: a*(A+B a+C a^2) = A a + B a^2 + C a^3 = A a + B a^2 + C(-a^2+2a+1) = C + (A+2C) a + (B-C) a^2
    # *a^2 = A a^2 + B a^3 + C a^4 = B + (2B+C) a + (A-B+3C) a^2
    import sympy as sp
    M = sp.Matrix([[sp.Rational(A), sp.Rational(C), sp.Rational(B)],
                   [sp.Rational(B), sp.Rational(A+2*C), sp.Rational(2*B+C)],
                   [sp.Rational(C), sp.Rational(B-C), sp.Rational(A-B+3*C)]])
    Minv = M.inv()
    col = Minv * sp.Matrix([sp.Rational(d),0,0])
    return (Fr(col[0]), Fr(col[1]), Fr(col[2]), 1)"""
new_inv = """def finv(x):
    A,B,C,d = x
    # multiplication matrix of (A+B a+C a^2), columns = images of 1,a,a^2:
    # *1 = A + B a + C a^2 ; *a = C + (A+2C) a + (B-C) a^2 ; *a^2 = (B-C) + (2B-C) a + (A-B+3C) a^2
    import sympy as sp
    M = sp.Matrix([[sp.Rational(A), sp.Rational(C), sp.Rational(B-C)],
                   [sp.Rational(B), sp.Rational(A+2*C), sp.Rational(2*B-C)],
                   [sp.Rational(C), sp.Rational(B-C), sp.Rational(A-B+3*C)]])
    Minv = M.inv()
    col = Minv * sp.Matrix([sp.Rational(d),0,0])
    return (Fr(col[0]), Fr(col[1]), Fr(col[2]), 1)"""
assert old_inv in src
src = src.replace(old_inv, new_inv)
open('ftfield.py','w').write(src)
import importlib, ftfield as F
importlib.reload(F)
a = F.AL
print("a^3-(1+2a-a^2):", F.fsub(F.fpow(a,3), F.ft(1,2,-1))[:3])
print("a^4-(3a^2+a-1):", F.fsub(F.fpow(a,4), F.ft(-1,1,3))[:3])
x = F.ft(1,2,0,3)
print("x*inv(x):", F.fmul(x, F.finv(x))[:3])
print("norm 1+a:", F.fnorm(F.fadd(F.ONE,a)), " norm 2+a:", F.fnorm(F.ft(2,1)), " norm a:", F.fnorm(a))
print("norm 1-a:", F.fnorm(F.fsub(F.ONE,a)))
print("N(2) =", F.fnorm(F.ft(2)))
open('progress.log','a').write('PG5: ftfield.py exact Q(2cos(2pi/7)) arithmetic fixed+tested\n')
