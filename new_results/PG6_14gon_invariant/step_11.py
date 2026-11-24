import pickle
base = '<project>/research/08_theory/'
pg5 = base+'PG5_alphabet_14gon/'
print("=== progress.log ===")
print(open(pg5+'progress.log').read())
for fn in ['gold14.pkl','E3n14_check.pkl','E3n14_spot.pkl']:
    print("="*40, fn)
    print(pickle.load(open(pg5+fn,'rb')))
