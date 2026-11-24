import shutil, os
SRC = '../CHN3_p5_locus'
files = ['point_engine.py','mzv_exact.py','kimilib.py','engine5.py','common.py',
         'C4.txt','Q7p.txt','subres1_A.txt','subres1_B.txt','R_a1a2a4.txt','R_factors.txt',
         'test_pt_w1112.pkl','forms10.pkl','linalg_tools.py']
for f in files:
    shutil.copy(os.path.join(SRC,f), f)
    print("copied", f, os.path.getsize(f))
# check engine5 imports
print("\n=== engine5.py ===")
print(open('engine5.py').read()[:1500])
