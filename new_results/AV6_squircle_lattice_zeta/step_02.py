import os, shutil
os.chdir("<project>/research/08_theory/AV6_squircle_lattice_zeta")
shutil.rmtree("research")
os.remove("progress.log")
with open("progress.log","w") as f:
    f.write("step_00: workspace init (cwd confirmed)\n")
print(os.getcwd(), os.listdir("."))
