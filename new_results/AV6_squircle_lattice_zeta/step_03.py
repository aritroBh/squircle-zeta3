import os
os.chdir("<project>/research/08_theory/AV6_squircle_lattice_zeta")
print(os.getcwd(), os.listdir("."))
with open("progress.log","w") as f:
    f.write("step_00: workspace init (cwd confirmed)\n")
print(os.listdir("."))
