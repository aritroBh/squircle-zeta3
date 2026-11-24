import pickle, sys, subprocess
sys.path.insert(0, '../DEC5D_decagon_order5_depth')
import dec5d_lib  # for QS5 unpickling
lb = pickle.load(open('../DEC5D_decagon_order5_depth/local_block_exact.pkl','rb'))
print(type(lb))
if isinstance(lb, dict):
    for k, v in lb.items():
        print(k, ':', str(v)[:400])
gold = pickle.load(open('../DEC5D_decagon_order5_depth/p5_gold.pkl','rb'))
print(gold)
