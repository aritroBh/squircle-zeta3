import os
wd = "<project>/research/08_theory/DEC2_decagon_exact_assembly"
files = sorted(os.listdir(wd))
key_files = ["verify_dec2.py","result.md","progress.log","geom_exact.pkl","E3n_sym_pq.pkl",
             "Jloc_dec2.pkl","Jtab_dec2.pkl","D_all_proved.pkl","MIX20_proved.pkl",
             "mons_dec2.pkl","S_assembly_dec2.pkl","final_table_dec2.pkl","esum_engine.py"]
for f in key_files:
    print(f, "OK" if f in files else "MISSING", os.path.getsize(os.path.join(wd,f)) if f in files else "")
print(open(os.path.join(wd,"progress.log")).read()[-700:])
