from o6util import load
import pickle
for name in ['atoms600.pkl','basis620.pkl','Jloc.pkl','rel_lll.pkl','P3_720.pkl']:
    d = load(name)
    print('='*20, name, type(d))
    if isinstance(d, dict):
        for k,v in d.items():
            print(' ', repr(k), '->', repr(v)[:160])
