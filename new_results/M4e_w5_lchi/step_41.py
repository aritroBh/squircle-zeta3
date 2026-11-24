import math
print("index factorization:")
n = 967458816000
import sympy as sp
print(sp.factorint(n))
# update nice_lchi.py -> nice_lchi2.py loading the saturated kernel
M4E = '<project>/research/08_theory/M4e_w5_lchi'
src = open(M4E+'/nice_lchi.py').read()
src = src.replace(
"""    K = [prim([F(c) for c in v]) for v in sym['kernel']]
    red = flint.fmpz_mat(K).lll(gram='exact')
    K = [[int(red[i, j]) for j in range(red.ncols())] for i in range(red.nrows())]""",
"""    ks = json.loads((HERE/'nice_ksat_s0.json').read_text())
    assert ks['names'] == names
    K = ks['K']""")
src = src.replace("lchi_search_{digits}.json", "lchi_sat_search_{digits}.json")
open(M4E+'/nice_lchi2.py','w').write(src)
print("written nice_lchi2.py (saturated kernel)")
