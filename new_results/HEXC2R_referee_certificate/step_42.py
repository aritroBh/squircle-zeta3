
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, ".")
sys.path.insert(0, "../COR4_M3alpha_classical")
import mpmath as mp
mp.mp.dps = 65
import reflib as R
from hpl_series import H_series
for w in [(0,), (1,), (-1,), (0,1), (0,-1), (1,-1), (0,0,1), (1,0,-1), (0,1,-1,1), (0,0,-1,1), (1,1,0,-1), (1,0)]:
    a = R.hpl_value(w, mp.mpf(1)/2); b = H_series(w, mp.mpf(1)/2)
    print(w, "err", mp.nstr(abs(a-b), 3), "| mine", mp.nstr(a, 8), "| theirs", mp.nstr(mp.mpc(b), 8))
