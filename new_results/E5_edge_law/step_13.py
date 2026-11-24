import os
os.makedirs('steps', exist_ok=True)
open('progress.log','a').write("E5 start. CWD is work dir. Reused E1: s01_edge_sym.py (radial-projection d-series), s03 beta table, D17 targets C_m^(5), C_m^(2,3).\n")
print(open('progress.log').read())
