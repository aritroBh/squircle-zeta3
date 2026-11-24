import os, subprocess, sys
print(os.getcwd())
os.makedirs("research/08_theory/AV6_squircle_lattice_zeta", exist_ok=True)
os.chdir("research/08_theory/AV6_squircle_lattice_zeta")
print(os.getcwd())
with open("progress.log","a") as f:
    f.write("step_00: workspace init\n")
print(os.listdir("."))
import mpmath
print("mpmath", mpmath.__version__)
print(sys.version)
