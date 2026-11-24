import pickle, os
base2 = "<project>/research/08_theory/Q1S_vertex_symbol/"
sols = pickle.load(open(base2+"sunit_li3.pkl","rb"))
print(type(sols), len(sols))
print(sols[0])
print(sols[1])
