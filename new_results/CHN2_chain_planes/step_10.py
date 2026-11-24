
import pickle, os
chn1 = '<project>/research/08_theory/CHN1_chain_bodies'
print(open('step_01.py').read())
print("="*60)
forms = pickle.load(open(f'{chn1}/forms.pkl','rb'))
print(type(forms), list(forms.keys())[:20] if hasattr(forms,'keys') else len(forms))
