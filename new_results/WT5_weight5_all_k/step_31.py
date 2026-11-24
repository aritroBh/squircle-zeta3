import shutil, os
base='<project>'
T=base+'/research/08_theory/'
W=T+'WT5_weight5_all_k/'
os.chdir(W)
copies = [
 (T+'WT4_weight4_all_k/', ['fastred.py','buildintegrand.py','genintegrand.py']+[f'binteg_{k}.pkl' for k in range(5,15)]),
 (T+'M5b_w6_closed_form/', ['w4_symbol.py','w5_symbol.py','w5_cobracket.py','w4_right_coproduct.py',
   'w5_search.py','w5_precision.py','g_eval_lib.py','pool.json','Lchi8_5.json','Lchi8_5_1800_planner.json',
   'omega5.json','gpl_w5pp.tsv','quad5_100.json','cobracket5pp.json','w5pp_closed_form.json','symbol5pp.json',
   'W5pp_1350.json','w5pp_4800.json','w5pp_vals_1450.pkl']),
]
for src, files in copies:
    for f in files:
        shutil.copy(src+f, W+f)
        print('copied', f)
with open(W+'progress.log','w') as f:
    f.write('step_01 setup: copied fastred/buildintegrand/genintegrand+binteg_5..14 (WT4), w4/w5 symbol+cobracket+precision+search+g_eval, pool/Lchi8/omega5/gpl_w5pp/quad5_100/cobracket5pp/w5pp_closed_form/W5pp_1350/w5pp_4800/w5pp_vals_1450 (M5b)\n')
print('files now:', sorted(os.listdir(W)))
