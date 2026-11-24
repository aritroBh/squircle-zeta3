import os
dst = '<project>/research/08_theory/M4i_w5_noL'
with open(dst+'/progress.log','a') as fh:
    fh.write('step_35(final): all checks done; verdict L droppable; L_relation 26 terms max|num|4815 lcm638976; R5_noL 69 terms max|num|2487169 lcm301916160 verified -2.1052657e-1596.\n')
print(open(dst+'/progress.log').read())
