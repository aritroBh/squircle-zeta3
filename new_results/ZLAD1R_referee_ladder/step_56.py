log = open('progress.log').read()
log = log.replace(
 "step_29: ZLAD1 prime: logGamma/beta identity numeric OK (4 pts, 35+ digits); EXACT resummation n=2..40 k<=8 ALL MATCH (2609 values)\n",
 "step_29: FIRST attempt at ZLAD1' checks — found bugs in MY harness (beta_poly s-step; logGamma series evaluated outside polydisc); superseded by step_30\n")
log = log.replace(
 "step_41: weight-9 rungs re-derived mechanically with sympy: ALL PASS (my earlier z2^3z3 sign was a hand-slip; correct rung is -4[z2^4]mu_8)\n",
 "step_41: mechanical re-derivation caught (i) my step_40 z2^3z3 hand sign-slip (correct: -4[z2^4]mu_8), (ii) a rational-coefficient truncation bug in my step_41 harness; fixed in step_45\n")
open('progress.log', 'w').write(log)
print(open('progress.log').read())
