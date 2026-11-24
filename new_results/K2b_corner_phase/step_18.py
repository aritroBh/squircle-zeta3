
import time
wd = 'research/08_theory/K2b_corner_phase/'
def log(msg):
    with open(wd+'progress.log','a') as f:
        f.write(time.strftime('%H:%M:%S ') + msg + '\n')
log("step_13: setup. Findings: (1) tables have ~65 sig digits; (2) free LS fits of kappa1 are biased (collinearity + image saddles e^{-(sqrt2-1)theta}); got 0.7776 vs true 0.8412 at pi/2 -> must derive analytically; (3) exact outer weight is e^{-x}/(1+e^{-x})^2 = sum_m (-1)^{m+1} m e^{-mx} -> dominant m=1 saddle + image towers; (4) plan: sympy saddle expansion for general alpha, match D13b at pi/2, then verify vs tables.")
print(open(wd+'progress.log').read())
