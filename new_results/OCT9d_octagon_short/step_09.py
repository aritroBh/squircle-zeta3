import pickle, sys
base = "<project>/research/08_theory/"
O9B = base+"OCT9b_octagon_debug"
s0 = pickle.load(open(O9B+'/sol_B0.pkl','rb'))
print(type(s0), list(s0.keys()))
for kk, vv in s0.items():
    print(kk, type(vv), str(vv)[:400])
