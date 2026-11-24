import pickle, sys
base = "<project>/research/08_theory/"
O9B = base+"OCT9b_octagon_debug"
sols2357 = pickle.load(open(O9B+'/sols2357.pkl','rb'))
print(type(sols2357))
print(str(sols2357)[:2000])
