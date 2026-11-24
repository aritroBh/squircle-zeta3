import os, subprocess
os.chdir('<project>')
print("CWD:", os.getcwd())
print("\n=== Work dir ===")
print(subprocess.run(['ls','-la','research/08_theory/FLD1_field_rule_proof'],capture_output=True,text=True).stdout)
print("\n=== G1_zeta3_any_angle ===")
print(subprocess.run(['ls','-la','research/10_frontier/G1_zeta3_any_angle'],capture_output=True,text=True).stdout)
print("\n=== W1 b3 dir ===")
print(subprocess.run(['ls','-la','research/10_frontier/W1_corner_family/B_taylor_pi4'],capture_output=True,text=True).stdout)