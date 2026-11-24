Q='<project>/research/08_theory/A3f_level8_relations/Q'
print(open(Q+'/check_conversion.py').read())
import pickle
xy,relations,cols,outputs=pickle.load(open(Q+'/kummer_aux.pickle','rb'))
print('kummer_aux: xy instances',len(xy),'relations',len(relations),'cols',len(cols),'outputs(eliminated H rows)',len(outputs))
print('nonzero outputs:',sum(bool(r) for r in outputs))
good=pickle.load(open(Q+'/kummer4.pickle','rb'))
print('kummer4 direct instances:',len(good))
