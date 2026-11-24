import pickle
cr = pickle.load(open('../PGEN1_polygon_alphabet/CR_presentations.pkl','rb'))
print(type(cr), list(cr.keys()) if isinstance(cr, dict) else len(cr))
for k in (cr if isinstance(cr, dict) else range(len(cr))):
    v = cr[k]
    print('edge', k, ':', str(v)[:400])
