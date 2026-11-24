
import pickle
chn1 = '<project>/research/08_theory/CHN1_chain_bodies'
forms = pickle.load(open(f'{chn1}/forms.pkl','rb'))
ks = sorted(forms.keys())
print(ks[:5], "...", ks[-5:])
v = forms[(2,0)]
print(type(v))
if isinstance(v, dict):
    items=list(v.items())[:5]
    print(items)
