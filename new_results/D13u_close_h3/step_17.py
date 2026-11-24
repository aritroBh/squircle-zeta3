import sys, math, time
sys.path.insert(0, '<project>/research/08_theory/D13s_fix_recertify/D13s_work/fixed')
import importlib
st = importlib.import_module('step_38')
t0=time.time()
r = st.certified_y38(math.sqrt(math.pi/(2*105.0)), math.sqrt(math.pi/(2*100.0)), J=24, nsub=4)
print(f"cell [100,105]: C={float(r['C'].upper()):.6f}  ({time.time()-t0:.1f}s)")
print(f"  -> estimated full step_18: ~{(time.time()-t0)*29/60:.1f} min")
