# Planner check of CX4's counterexample to WT2b's displayed Lemma B (key (2,1), k = 4).
import sympy as sp
I, r2 = sp.I, sp.sqrt(2)
t1, t2 = (-1-I)/2, (-1+I)/2
h31, h22, h21_actual = -5*r2*I/4, 15*r2*I/16, -3*r2/4
adj = 2*t2*(h22 - sp.conjugate(h22))
printed = (3*t1*h31 + 2*t2*h22 + adj)/5      # [i>j] term plus the adjacent term: the diagonal feed counted twice
corrected = (3*t1*h31 + adj)/5               # [i>=j+2] term only
print("printed  :", sp.simplify(printed), " actual - printed =", sp.simplify(h21_actual - printed))
print("corrected:", sp.simplify(corrected), " actual - corrected =", sp.simplify(h21_actual - corrected))
