import sys, json, os
sys.path.insert(0,'<project>/research/08_theory/M5b_w6_closed_form')
os.chdir('<project>/research/08_theory/M5b_w6_closed_form')
import m5b_precision as MP
try:
    d = MP.run('w6', 600)
    print('OK', d['radius'], d['tail_bound'])
except AssertionError as e:
    import traceback; traceback.print_exc()
