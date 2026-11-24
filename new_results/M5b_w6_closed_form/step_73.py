import sys
sys.path.insert(0,'<project>/research/08_theory/M5b_w6_closed_form')
import os
os.chdir('<project>/research/08_theory/M5b_w6_closed_form')
from w4_symbol import Q, R, A, B, g_symbol, add, scale
from m5b_common import ZERO, ONE
import w5pp_symbol as WS
from w6_solve import n6_terms, mpl_word, parse_qs2
for m,n,x,y,c in n6_terms():
    print(m,n,'x=',x,'y=',y,'c=',c)
    w = mpl_word(m,n,x,y)
    print('  word:', w)
    try:
        s = g_symbol(ZERO, w, ONE)
        print('  symbol ok, terms', len(s))
    except Exception as e:
        print('  FAIL:', e)
