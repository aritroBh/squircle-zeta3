base = '<project>/research/08_theory'
c = base + '/ZQ5c_zagier_w5_capped'
print('===== step_19_d13cap_vals.py'); print(open(c+'/step_19_d13cap_vals.py').read())
import pickle
raw = pickle.load(open(c+'/out/d13/vals_cap_2000.pkl','rb'))
print('vals keys sample:', list(raw.keys())[:12], '... total', len(raw))
