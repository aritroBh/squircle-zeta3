import sys, pickle
sys.path.insert(0, '../DEC5D_decagon_order5_depth')
import dec5d_lib as D
bs = pickle.load(open('../DEC5D_decagon_order5_depth/dec_block_symbols.pkl','rb'))
tot = bs['tot']
items = sorted(tot.items(), key=lambda kv: kv[0])
print("total words in S(W5/L):", len(items))
for w, c in items[:8]:
    print(w, c)
print('...')
sub = bs['sub_sym']
for k in sub:
    print(k, len(sub[k]))
