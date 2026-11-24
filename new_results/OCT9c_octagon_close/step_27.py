import pickle
z = pickle.load(open('../OCT9_generic_octagon_closed/Ze_check.pkl','rb'))
print(type(z))
import pprint
if isinstance(z, dict):
    for k,v in z.items(): print(k, ':', str(v)[:400])
else:
    print(z)
