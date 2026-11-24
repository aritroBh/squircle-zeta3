
import pickle
br = pickle.load(open('basis_red81.pkl','rb'))
print(type(br), list(br.keys()) if isinstance(br, dict) else len(br))
br2 = br['basis_red']
print(len(br2), type(br2[0]))
