import os, pickle
# confirm files present
need = ["q1s_exact.py","q1s_symbol.py","q3_engine.py","q3_eval.py","q3_stage2.py","q3_symbol.py","sunit_li3.pkl","q3_solve.pkl","q3_symbols.pkl","q3_accs.pkl"]
print({f: os.path.exists(f) for f in need})
sol = pickle.load(open('q3_solve.pkl','rb'))
for nm in sol:
    print(nm, "li3:", len(sol[nm]['li3']), "dec:", len(sol[nm]['dec']))
print(sol['V3p']['dec'])
sols = pickle.load(open('sunit_li3.pkl','rb'))
print(sols[145]['z'], sols[145]['vec_z'], sols[145]['vec_1z'])
