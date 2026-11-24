
import sys; sys.path.insert(0,'.')
import json, time
import mpmath as mp
import w5lib as W
import pipe
D=17
C = W.letters_info(D)
J, cands = pipe.load_cands(D)
t0=time.time()
vals = pipe.phase_numerics(D, C, cands, P=1150)
print(f"numerics {time.time()-t0:.0f}s")
mp.mp.dps=60
print("L(chi17,5) =", mp.nstr(vals["L"], 30))
print("T = L/sqrt17 =", mp.nstr(vals["T"], 30))
