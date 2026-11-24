import sys
sys.path.insert(0, "<project>/research/10_frontier/G1_zeta3_any_angle")
sys.path.insert(0, "<project>/research/08_theory/COR4_M3alpha_classical")
import mpmath as mp
from hpl_series import H_series
import hpl_eval as he

mp.mp.dps = 60
x = mp.mpf(1)/2
for w in [(0,1),(0,-1),(1,),(  -1,),(0,0,1),(0,1,-1),(1,0,1),(0,0,0,1),(0,1,0,-1),(1,0,-1,1),(0,0,1,-1)]:
    try:
        v_he = he.H(w, x)
    except Exception as ex:
        v_he = f"ERR {ex}"
    v_hs = H_series(w, x)
    print(w, "hpl_eval:", mp.nstr(v_he, 25) if not isinstance(v_he, str) else v_he)
    print(w, "H_series:", mp.nstr(v_hs, 25))
