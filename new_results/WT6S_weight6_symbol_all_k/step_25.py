
import os, sys
WD = '<project>/research/08_theory/WT6S_weight6_symbol_all_k'
os.chdir(WD)
sys.dont_write_bytecode = True
if WD not in sys.path:
    sys.path.insert(0, WD)
import qsqrt2  # import local first so it's cached
from w4_symbol import (Q, R, A, B, ZERO, ONE, g_symbol, log_product, LT, MT,
                       shuffle, log_symbol, add, scale, concat, clean, ban_float, self_check)
print("import OK")
print("A =", A, " B =", B, " R =", R)
# quick check: g_symbol of a single letter
s = g_symbol(ONE, (ZERO,), A)
print("G_1(0;a) symbol:", s)
# LT structure
print("LT:", LT)
print("MT:", MT)
