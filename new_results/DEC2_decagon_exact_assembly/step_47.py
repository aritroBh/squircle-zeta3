import os, subprocess
# search for where Jtab_t.pkl is built
hits = subprocess.run(['grep','-rl','Jtab_t', '<project>/research/08_theory/K1d_engine_fixed'], capture_output=True, text=True).stdout
print(hits)
hits2 = subprocess.run(['grep','-rln','esum', '<project>/research/08_theory/K1d_engine_fixed'], capture_output=True, text=True).stdout
print(hits2)
