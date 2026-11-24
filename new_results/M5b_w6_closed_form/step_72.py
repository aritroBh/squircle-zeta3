import sys, itertools
sys.path.insert(0,'<project>/research/08_theory/M5b_w6_closed_form')
import os
os.chdir('<project>/research/08_theory/M5b_w6_closed_form')
from w4_symbol import Q, R, A, B, g_symbol, li_symbol, li_symbol as ls
from m5b_common import ZERO, ONE, li32_word, W4_SEVEN
from fractions import Fraction as F
# test each li32 atom
for x,y in [(B,ONE),(R/2,ONE),(-B,ONE),(B,-ONE)]:
    try:
        s=g_symbol(ZERO, li32_word(x,y), ONE)
        print('li32 ok', x, y)
    except Exception as e:
        print('li32 FAIL', x, y, e)
# test li6 args
for z in W4_SEVEN + [Q(F(1,2))-R/2, Q(F(1,2))-R/4, Q(F(1,2))+R/4]:
    try:
        s=li_symbol(6,z)
        print('li6 ok', z)
    except Exception as e:
        print('li6 FAIL', z, e)
