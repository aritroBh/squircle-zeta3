
import sympy as sp, pickle
WD = '<project>/research/08_theory/K1f_hexagon_order5'
K1D = '<project>/research/08_theory/K1d_engine_fixed'
e3 = pickle.load(open(f'{K1D}/E3_fixed.pkl','rb'))
print("keys:", e3.keys())
print("nl count:", len(e3['nl']), "loc count:", len(e3['loc']))
print("sample nl:", e3['nl'][:3])
print("sample loc:", e3['loc'][:2])
