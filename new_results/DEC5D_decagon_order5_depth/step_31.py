
import pickle
base = "<project>/research/08_theory/"
ft = pickle.load(open(base+"DEC1_decagon/final_table_pq.pkl","rb"))
print(type(ft))
print(ft)
