import os
# inventory of key artifacts for the final report
arts = ['result.md','progress.log','mons_simple.pkl','Dconsts.pkl','mixed_engine.py','S_assembly.pkl',
        'reductions.pkl','basis745.pkl','indep.pkl','S_indep.pkl','final_tables.pkl','S_final.pkl',
        'table_Stot.pkl','table_Smix.pkl','table_diffC4.pkl','subtables.pkl','P3_exact.pkl','num745.pkl']
for a in arts:
    print(f"{a:20s}", os.path.getsize(a) if os.path.exists(a) else "MISSING")
print()
print(open('progress.log').read())
