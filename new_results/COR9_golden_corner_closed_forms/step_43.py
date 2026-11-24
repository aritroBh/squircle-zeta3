import sympy as sp
WD = '<project>/research/08_theory/COR9_golden_corner_closed_forms/'
names = ["zeta3", "L1=Li3(1/phi)", "L2=Li3(-1/phi)", "L3=Li3(1/phi^3)", "L4=Li3(-1/phi^3)",
         "L5=Li3(1/(2phi))", "L6=Li3(2/phi^2)", "L7=Li3(1/(2phi^2))",
         "pi2*L", "pi2*lam", "Li2(1/phi^3)*L", "Li2(1/phi^3)*lam", "L^3", "L^2*lam", "L*lam^2", "lam^3"]
ra = [0, 6,-6,-8, 8, 0, 0, 0, -1, 0, 0, 0, 9, 0, 0, 0]
rb = [3,-18,-12, 4,-4, 0, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0]
rc = [18,-3, 18, 6,-22, 0,-12,-12, -2, 1, 0, 0, 35,-48, 12, 2]
rd = [21,-15, 0, 2,-18, 0,-12,-12, -3, 1, 0, 0, 42,-48, 12, 2]
R = sp.Matrix([ra, rb, rc, rd])
print("rank of relation matrix:", R.rank())
pivots = R.rref()[1]
print("pivot columns:", pivots, "-> eliminate these elements:", [names[i] for i in pivots])
free = [i for i in range(16) if i not in pivots]
print("independent basis (free columns):", [names[i] for i in free])
