import pickle
names=pickle.load(open('q3_names.pkl','rb'))['names']
sel=pickle.load(open('q3_sel.pkl','rb'))['sel']
solmin=pickle.load(open('q3_solve_min.pkl','rb'))['V3p']
print(type(names[0]), names[0])
print(type(sel[0]), sel[0])
print(type(solmin[0]), solmin[0])
print([type(v) for k,v in names[:3]], [type(v) for k,v in sel[:3]])
