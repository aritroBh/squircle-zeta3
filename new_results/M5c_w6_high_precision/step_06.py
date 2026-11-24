import os, shutil
base='<project>/research/08_theory/M5b_w6_closed_form'
work='<project>/research/08_theory/M5c_w6_high_precision'
scripts=['m5b_precision.py','w5pp_symbol.py','g_eval_lib.py','g_eval5pp.py','w5pp_cobracket.py','m5b_common.py',
'w5pp_solve.py','w5pp_lll.py','w5pp_assemble.py','w5pp_verify.py','w6_solve.py','w6_values.py','w6_lll.py',
'w6_lll_final.py','w6_assemble.py','w6_verify.py','m5b_indep_check.py',
'w5_symbol.py','w5_cobracket.py','w5_search.py','w5_precision.py','nice_search.py','w4_symbol.py','w4_right_coproduct.py']
for s in scripts:
    shutil.copy(os.path.join(base,s), os.path.join(work,s))
# also copy data jsons that are inputs (read-only use)
data=['W6_5400.json','Lchi8_5_2500.json','w6_symbol_solve.json','n6_solution.json','R5_min3.json','pool.json']
for s in data:
    shutil.copy(os.path.join(base,s), os.path.join(work,s))
print(sorted(os.listdir(work)))
