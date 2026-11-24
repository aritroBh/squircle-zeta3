import json
base = '<project>/research/08_theory'
c = base + '/ZQ5c_zagier_w5_capped'
gl = json.load(open(c+'/out/d13/genuine_lattice_cap.json'))
G = [[int(x) for x in row] for row in gl['G']]
print('rows:', len(G), 'cols:', len(G[0]))
names = gl['atom_names']
print('n names:', len(names))
print(names[:40])
print(names[40:])
for i,row in enumerate(G):
    print('row',i,'T_coef=',row[0], 'z5=',row[1], 'nnz=',sum(1 for x in row if x))
