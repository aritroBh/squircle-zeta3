import shutil, os
cwd = os.getcwd()
src = '../CHN2_chain_planes'
for f in ['engine5.py','kimilib.py','mzv_exact.py','mzv_num.py','mzv_vals.py','symengine.py','engine.py',
          'forms10.pkl','forms8_W8.pkl','fft_check.pkl','numcheck_pt1.pkl','numcheck_pt2.pkl']:
    shutil.copy(os.path.join(src,f), os.path.join(cwd,f))
# also from CHN1
src1 = '../CHN1_chain_bodies'
for f in ['w12_planes.json','w12_planes.py']:
    shutil.copy(os.path.join(src1,f), os.path.join(cwd,f))
print("copied. my dir now:")
for f in sorted(os.listdir(cwd)): print("  ", f)
print("\n--- CHN2 progress.log ---")
print(open(os.path.join(src,'progress.log')).read())
print("\n--- planner dir ---")
for f in sorted(os.listdir(os.path.join(src,'planner'))): print("  ", f)
