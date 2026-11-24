
import os, json
WD = '<project>/research/08_theory/WT6S_weight6_symbol_all_k'
os.chdir(WD)
Ev = json.load(open('e5m_values.json'))['values']
for m in range(1,6):
    print(f"E5,{m} = {Ev[f'E5,{m}'][:55]}")
ns = json.load(open('n6_symbols_pairings.json'))
print("\nsymbol words:", ns['symbol_words'])
print("\nN6 polynomial (depth-two of W6(k)):")
np_ = json.load(open('n_w6k_polynomial.json'))
for name, v in np_.items():
    print(f"  [{v['rat']} + ({v['sqrt2']})*sqrt2] * {name}")
