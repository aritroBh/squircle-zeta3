import pickle
d='../OCT9b_octagon_debug/'
lists = pickle.load(open(d+'Sunit_lists.pkl','rb'))
print(type(lists), list(lists.keys()))
G = lists['Gargs']; D = lists['Dargs']
print(len(G), len(D))
print('Gargs:', G)
print('Dargs:', D)
