import os
os.chdir('<project>/research/08_theory/ZQ5b_zagier_w5_larger_S')
with open('progress.log','a') as f:
    f.write('step_03a_d17_crosscheck: my factor_sunit reproduces ZQ5 d17 vx/v1m on 40/40 sampled candidates\n')
    f.write('step_03b_d17_verify: planner d17 relation independently verified: symbol ZERO, cobracket ZERO (full/cut, total/Li5), fresh resid -5.96e-1319 @1300d; 94 Li5+34 prod+z5; T_coef=-17^5; L=1.02763403814516292324626545432\n')
    f.write('step_04a_d17_setup: relation expanded to atom space (matches JSON), 43-dim LLL-reduced zero lattice saved\n')
    f.write('step_04b_d17_vals: values @1200d; relation resid -6.7e-1219; worst zero-row resid 1.375e-1220\n')
    f.write('step_04e_d17_minA: max numerator 232217 -> 89612 (103 Li5, 137 terms)\n')
    f.write('step_04f_d17_minB: term-min form 74 Li5 / 1630051 / 109 terms (frontier)\n')
    f.write('step_05_d17_final: d17_final.json (maxnum-min): symbol ZERO, cobracket ZERO, fresh resid -3.74e-1319 @1300d; stats 94/232217, 103/89612, 74/1630051\n')
print(open('progress.log').read())
