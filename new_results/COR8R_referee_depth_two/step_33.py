import cor8r_sym as S
import sympy as sp
cc = S.cc
A = 1/(1+cc); B = 2/(1+cc)
# try bad words of increasing length
tests = [
 (A, sp.Integer(0)),                    # L2 trailing zero
 (sp.Integer(0), A, sp.Integer(0)),     # L3 trailing zero
 (sp.Integer(0), sp.Integer(0), A),     # L3 base (leading zeros) - should use base case
 (A, sp.Integer(0), sp.Integer(0)),     # L3 interior repeat
 (sp.Integer(0), sp.Integer(0), A, B),  # L4 leading repeat
 (A, sp.Integer(0), sp.Integer(0), B),  # L4 interior repeat
]
for w in tests:
    try:
        s = S.symG(w)
        print("OK  ", w, "->", len(s), "terms")
    except Exception as e:
        print("FAIL", w, "->", type(e).__name__, e)
