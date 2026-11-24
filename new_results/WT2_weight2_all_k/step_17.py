import sys, time
sys.dont_write_bytecode = True
sys.path.insert(0, "<project>/paper/squircle_note/scripts")
import sympy as sp
from mk_integrands import t_integrand, v_integrand

for k in [6,7,8]:
    t0=time.time(); tk = t_integrand(k); t1=time.time()
    print(f"k={k}: t_integrand {t1-t0:.2f}s, terms={len(tk)}", flush=True)
    t0=time.time(); vk = v_integrand(k); t1=time.time()
    print(f"      v_integrand {t1-t0:.2f}s", flush=True)
