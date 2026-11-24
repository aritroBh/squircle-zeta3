import pickle
d='../OCT9b_octagon_debug/'
s0 = pickle.load(open(d+'sol_B0.pkl','rb'))
meta = s0['meta']
print(len(meta['sols']), len(meta['negw']), len(meta['logprods']))
print(len(s0['solB0']), 375+len(meta['negw'])*4+20)
print(open(d+'step_83.py').read()[:6000])
