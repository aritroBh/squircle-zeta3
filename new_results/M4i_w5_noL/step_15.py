import shutil, os, time
base = '<project>'
m4f = base + '/research/08_theory/M4f_w5_extras'
dst = base + '/research/08_theory/M4i_w5_noL'
os.makedirs(dst, exist_ok=True)
files = ['m4f_common.py','nice_search.py','w5_symbol.py','w5_search.py','w5_cobracket.py',
         'w5_precision.py','w5_verify.py','overcomplete.py','overcomplete_1200.json',
         'overcomplete_1500.json','subset_exact.py','subset_exact.pkl','drop_L.py',
         'optimize6.py','optimize7.py','export_forms.py','verify_json.py',
         # inputs
         'pool.json','nice_ksat_s0.json','nice_symbol.json','nice_atoms.pkl','cobracket.json',
         'Lchi8_5.json','W5_5400.json','R5_best.json','R5_pareto.json','R5_min3.json',
         'state_opt.pkl','state_opt2.pkl','pareto.pkl','state_min3.pkl',
         'vals_cache_1200.pkl','vals_cache_1500.pkl','height_lll.py']
for f in files:
    shutil.copy2(os.path.join(m4f,f), os.path.join(dst,f))
# record as step_01
open(os.path.join(dst,'step_01.py'),'w').write(
    '"""step_01: copy (never edit in place) M4f scripts+data into M4i_w5_noL.\n'
    'Files: '+', '.join(files)+'\n"""\n')
with open(os.path.join(dst,'progress.log'),'a') as fh:
    fh.write('step_01: copied %d M4f scripts+data files into M4i_w5_noL (M4f untouched).\n' % len(files))
print('copied', len(files), 'files to', dst)
print(sorted(os.listdir(dst)))
