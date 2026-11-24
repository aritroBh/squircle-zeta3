import os, time
t5 = '<project>/research/08_theory/TOW5_tower_depth10'
print(open(t5+'/progress.log').read())
print('======= verify_tow5.log =======')
print(open(t5+'/verify_tow5.log').read())
files = ['tow5lib.py','prove_S5.py','N5_proved.pkl','R5_raw.pkl','tau910.pkl','L910_candidate.pkl','a9_closed.pkl','a10_closed.pkl','forced_F910.pkl','forced910.py','data_rho_D10_d7_40.pkl','data_towers_D10_K24.pkl','data_towers_D10_K30_deg18.pkl','verify_tow5.py','prove_adj910_fast.py','prove_E3b910.py','check_full910.py','verify_W7.py','taskE_defect.py','L_simplified.pkl','L_parsed.pkl','forced_F78.pkl']
for f in files:
    p = os.path.join(t5,f)
    print(f, time.strftime('%m-%d %H:%M', time.localtime(os.path.getmtime(p))))
