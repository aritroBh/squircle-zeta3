import os
K1d="research/08_theory/K1d_engine_fixed"
print(os.path.exists(K1d))
if os.path.exists(K1d):
    for f in sorted(os.listdir(K1d)):
        print(f, os.path.getsize(os.path.join(K1d,f)))
