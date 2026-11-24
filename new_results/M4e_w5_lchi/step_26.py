import os, shutil, sys
print("cwd:", os.getcwd())
print("python:", sys.executable)
M4E = '<project>/research/08_theory/M4e_w5_lchi'
M4D = '<project>/research/08_theory/M4d_w5_nice'
for f in ['w5_symbol.py','w5_search.py','w5_precision.py','nice_search.py',
          'W5_5400.json','cobracket.json','nice_symbol.json','nice_atoms.pkl','nice_ksat.json']:
    shutil.copy(os.path.join(M4D,f), os.path.join(M4E,f))
print(sorted(os.listdir(M4E)))
