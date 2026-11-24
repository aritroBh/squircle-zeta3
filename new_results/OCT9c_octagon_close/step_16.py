import pickle, sympy as sp
d='../OCT9b_octagon_debug/'
def L(name): return pickle.load(open(d+name,'rb'))
s123 = L('sol_B123.pkl')
s0 = L('sol_B0.pkl')
for k,v in s123.items():
    print('==== sol B', k, type(v))
    if isinstance(v, dict):
        for kk, vv in v.items():
            print('   ', kk, type(vv), (str(vv)[:300] if not isinstance(vv,(list,tuple)) else f'len {len(vv)}'))
    elif isinstance(v,(list,tuple)):
        print('len', len(v)); print(v[:5])
print()
print('==== sol_B0 meta:', s0['meta'])
solB0 = s0['solB0']
print(type(solB0), len(solB0) if hasattr(solB0,'__len__') else '')
