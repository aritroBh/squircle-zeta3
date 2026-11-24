import pickle
d = pickle.load(open('../OCT6_octagon_closed/final_tables.pkl','rb'))
tot = d['tot']
for w, nm, cf in tot:
    print(w, nm, '=', cf)
