M4E = '<project>/research/08_theory/M4e_w5_lchi'
# quick sanity: list produced files
import os
fs = ['Lchi8_5.json','nice_sat_s0.py','nice_ksat_s0.json','nice_lchi.py','nice_lchi2.py',
      'lchi_diag.py','lchi_diag_800.json','lchi_second.py','lchi_testB.py','lchi_testB_1200.json',
      'lchi_min.py','lchi_min_1200.json','lchi_final.py','lchi_final_1200.json','lchi_final_1500.json',
      'verify_nice.py','closed_form_nice.md','assemble_m4_nice.py','n5_factor7.py','lchi_valrel.py',
      'lchi_sat_search_1200.json','lchi_sat_search_1500.json','progress.log']
for f in fs:
    print(f, os.path.exists(os.path.join(M4E,f)), os.path.getsize(os.path.join(M4E,f)) if os.path.exists(os.path.join(M4E,f)) else '-')
