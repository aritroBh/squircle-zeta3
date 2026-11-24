import pickle, sys
base = '<project>/research/08_theory/'
work = base+'TOW3R_review/'
ef = pickle.load(open(work+'tow3_copy/e_fits.pkl','rb'))
print('e_fits.pkl keys:', list(ef.keys())[:20] if hasattr(ef,'keys') else type(ef))
e61 = pickle.load(open(work+'tow3_copy/e61_closed.pkl','rb'))
print('e61_closed.pkl:', type(e61), list(e61.keys()) if hasattr(e61,'keys') else e61)
for k,vv in (ef.items() if hasattr(ef,'items') else []):
    print(k, ':', str(vv)[:400])
