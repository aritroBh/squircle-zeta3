import sympy as sp
from sympy.polys.matrices import DomainMatrix

# Test DomainMatrix
sqrt2 = sp.sqrt(2)
Qsqrt2 = sp.QQ.algebraic_field(sqrt2)

# Simple test
test_data = [[1 + sqrt2, 2], [3, 1 - sqrt2]]
print("Test data:", test_data)

try:
    dm = DomainMatrix.from_list(test_data, Qsqrt2)
    print("DomainMatrix created:", dm)
    print("Shape:", dm.shape)
    print("Rank:", dm.rank())
except Exception as e:
    print(f"Error: {e}")
    # Try different syntax
    try:
        rows = len(test_data)
        cols = len(test_data[0])
        dm = DomainMatrix.from_list_sympy(rows, cols, test_data, Qsqrt2)
        print("Alternative worked")
    except Exception as e2:
        print(f"Also failed: {e2}")