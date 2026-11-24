
from flint import ctx, arb, acb
ctx.dps = 40
from step_07_partb import *
print("S_of(400) =", S_of(400).str(10))
print("exp(S) =", S_of(400).exp().str(10))
print("H_of(400) =", H_of(400).str(10))
print("momentsum =", momentsum_of(400).str(10))
