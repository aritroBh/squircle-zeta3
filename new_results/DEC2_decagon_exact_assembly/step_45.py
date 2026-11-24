import os
base = "<project>/research/08_theory/K1d_engine_fixed"
for f in sorted(os.listdir(base)):
    if f.endswith('.py'):
        txt = open(os.path.join(base,f)).read()
        if 'Jtab' in txt and ('esum' in txt or 'pickle.dump' in txt):
            print("=====", f, "=====")
            print(txt[:4000])
            break
