import time, sys
sys.path.insert(0,'.')
import point_engine as PE, mzv_exact as M
from fractions import Fraction as F

# Time the exact Fraction engine at W11 at a small integer point
t0=time.time()
vw = PE.volreg_words_point((1,2,3,4,5), 11)
t1=time.time()
print(f"volreg_words_point W11 forward: {t1-t0:.2f}s, words at W11: {len(vw[11])}")
vwr = PE.volreg_words_point((5,4,3,2,1), 11)
t2=time.time()
print(f"volreg_words_point W11 reversed: {t2-t1:.2f}s")
# build table W11
ta=time.time()
T11=PE.table(11)
print(f"table(11) build: {time.time()-ta:.2f}s, nullity {T11.nullity}")
