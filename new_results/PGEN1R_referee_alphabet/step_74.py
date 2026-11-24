import pickle
T = pickle.load(open('../PGEN1_polygon_alphabet/field_tables.pkl','rb'))
print(T.keys())
T3 = T['T3']
for k in sorted(T3):
    print(k, '->', T3[k])
