
import time, pickle
from fractions import Fraction as F
import mzv_exact as M

t=time.time(); T11=M.DSTable(11); print(f"DSTable(11) pure: {time.time()-t:.1f}s nullity={T11.nullity}")
t=time.time(); T11f=M.DSTableFlint(11); print(f"DSTableFlint(11): {time.time()-t:.1f}s")
t=time.time(); T12f=M.DSTableFlint(12); print(f"DSTableFlint(12): {time.time()-t:.1f}s nullity={T12f.nullity} words={T12f.n}")
