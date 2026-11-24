cwd = "<project>/research/08_theory/G2_M4alpha"
hdr = '"""step_05.py -- EXACT symbolic formula for mu_{4,n} (k=4), extending mu3_symbolic.py to weight-4 zeta content.\n\n'
hdr += 'mu_{4,n} = -(1/(n(2n-1))) sum_{(j,i,l)} e_{j,i,l}(n) (n)_j/(2n)_j D_{i,l}(n+j, n),\n'
hdr += 'as an exact polynomial in h1..h4 = H^(r)_{n-1}, g1..g4 = H^(r)_{2n-1}, z2, z3, z4 with Q(n) coefficients.\n'
hdr += 'New at k=4: the third polygamma gives 6 z4 - 6 H4; monomials z2^2 and z4 occur (both prop. to pi^4).\n'
hdr += 'Checks: (1) formula vs flat_exact.py direct (n<=30); (2) patterns for the z4, z2^2, z3 parts;\n'
hdr += 'Writes mu4_formula.pkl and flat_exact_N400_K4_formula.json (G1 hpl_fit format).\n"""\n'
body = open(cwd + "/step_05.py").read()
# strip the old broken docstring (first triple-quoted block)
i = body.find('"""')
j = body.find('"""', i + 3)
rest = body[j + 3:]
code = hdr + rest
with open(cwd + "/step_05.py", "w") as f:
    f.write(code)
import ast
ast.parse(code)
print("syntax ok, length", len(code))
