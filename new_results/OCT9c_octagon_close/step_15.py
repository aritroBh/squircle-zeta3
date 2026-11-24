import pickle, sympy as sp
d='../OCT9b_octagon_debug/'
def L(name): return pickle.load(open(d+name,'rb'))
s123 = L('sol_B123.pkl')
s0 = L('sol_B0.pkl')
print(type(s123), list(s123.keys()) if isinstance(s123,dict) else '')
print(type(s0), list(s0.keys()) if isinstance(s0,dict) else '')
