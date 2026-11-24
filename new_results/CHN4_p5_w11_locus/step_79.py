import pickle, sympy as sp
forms=pickle.load(open('forms10.pkl','rb'))
print("type:", type(forms))
if isinstance(forms, dict):
    print("keys:", list(forms.keys()))
elif isinstance(forms,(list,tuple)):
    print("len:", len(forms))
    for i,x in enumerate(forms):
        print(i, type(x), (list(x.keys())[:5] if isinstance(x,dict) else x))
