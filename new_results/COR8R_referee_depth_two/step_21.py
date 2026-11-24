import cor8r_sym as S
import sympy as sp
cc = S.cc
# factorization tests
print(S.factor_Qc(2))                    # [2]
print(S.factor_Qc(sp.Rational(1,2)))     # -[2]
print(S.factor_Qc(-cc))                  # [-1]+[c]
print(S.factor_Qc((1-cc)/(2*cc)))        # [1-c]-[2]-[c]
print(S.factor_Qc(2/(1+cc)))             # [2]-[1+c]
print(S.factor_Qc((cc-1)/(1+cc)))        # [1-c]+[-1]-[1+c]
# base case G(A;1) = log(1-1/A)
A = 2/(1-cc)
print("S(G(A;1)):", S.symG((A,)))
print("S(G(0,A;1)):", S.symG((sp.Integer(0),A)))
print("S(G(0,0,A;1)):", S.symG((sp.Integer(0),sp.Integer(0),A)))
