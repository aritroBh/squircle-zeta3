import os
print("CWD:", os.getcwd())
print("\n--- E5d_nonlocal_pslq ---")
for f in sorted(os.listdir('../E5d_nonlocal_pslq')):
    print('  ', f)
print("\n--- E5c_nonlocal_z5 ---")
for f in sorted(os.listdir('../E5c_nonlocal_z5')):
    print('  ', f)
print("\n--- CWD ---")
for f in sorted(os.listdir('.')):
    print('  ', f)
