import shutil, os
SRC = 'research/08_theory/K4c_3d_edge_law/'

# Copy driver, pslq lib, engine, quad, and data files
for f in ['k4c_driver.py','pslq_lib.py','k4_engine.py','k4_quad.py',
          'modelA_h32.json','modelA_h48.json','modelA_h64.json','modelA_h96.json',
          'modelB_h32.json','modelB_h48.json','modelB_h96.json',
          'K1V_A.npy','K1V_B.npy']:
    shutil.copy(os.path.join(SRC,f), f)
    print("copied", f)

# Read existing step files
print("\n=== existing step_01.py ===")
print(open('step_01.py').read())
print("=== existing step_02.py ===")
print(open('step_02.py').read())
