import sys; sys.path.insert(0,'.')
import numpy as np, importlib
import modengine as ME, point_engine as PE, engine5, kimilib
importlib.reload(ME)
p=2147483647
a=(1,1,1,1,1)
# engine5 symbolic reference for Est, Inn, Gg
eng = engine5.build_engine(5)
for key in ['Est','Inn','Gg']:
    print(f"--- engine5 {key} ---")
    for j in range(2,6):
        d=eng[key][j]
        if d: print(f" {key}[{j}]:", {w: kimilib.peval5(pp,a) for w,pp in d.items()})
