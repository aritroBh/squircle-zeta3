import pickle, sys
base = "<project>/research/08_theory/"
O9B = base+"OCT9b_octagon_debug"
s0 = pickle.load(open(O9B+'/sol_B0.pkl','rb'))
meta = s0['meta']
print("meta keys:", list(meta.keys()))
print("n sols:", len(meta['sols']), "n negw:", len(meta['negw']), "n logprods:", len(meta['logprods']))
print("sols:", meta['sols'])
print("negw:", meta['negw'])
print("logprods:", meta['logprods'])
print("solB0 length:", len(s0['solB0']))
lists = pickle.load(open(O9B+'/Sunit_lists.pkl','rb'))
print("lists keys:", list(lists.keys()))
print("Gargs (%d):"%len(lists['Gargs']), lists['Gargs'])
print("Dargs (%d):"%len(lists['Dargs']), lists['Dargs'])
