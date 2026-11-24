from mpmath import mp, mpf
mp.dps = 60
import sys, time, json, os; sys.path.insert(0,'.')
import importlib, squircle; importlib.reload(squircle)
from squircle import bn_der
t0=time.time(); v=bn_der(1,1); print("b1' :", mp.nstr(v,40), f"{time.time()-t0:.1f}s")
t0=time.time(); v=bn_der(3,3); print("b3''':", mp.nstr(v,40), f"{time.time()-t0:.1f}s")
