import shutil, os, sys
base='<project>/research/08_theory/A3f_level8_relations'
work='<project>/research/08_theory/A3f2_level8_w3w4'
# copy Q dir
if not os.path.isdir(work+'/Q'):
    shutil.copytree(base+'/Q', work+'/Q')
for f in ['audit_weight1.pickle','audit_weight2.pickle','audit_weight3.pickle',
          'weight2_new_relation.pickle','weight2_new_relation.json','weight2_fiveterm.pickle',
          'weight3_with_lifts.pickle','weight3_onevar.pickle']:
    if os.path.exists(base+'/'+f) and not os.path.exists(work+'/'+f):
        shutil.copy(base+'/'+f, work+'/'+f)
print(sorted(os.listdir(work)))
# check qsqrt2 location
print(os.path.exists('<project>/paper/squircle_note/scripts/qsqrt2.py'))
print(os.path.exists('<project>/research/08_theory/A3_W4_proof/w4_symbol.py'))
