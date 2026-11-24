
import os
print("CWD:", os.getcwd())
print(os.listdir("."))
print("---- read-only DEC5D:")
rd = "<project>/research/08_theory/DEC5D_decagon_order5_depth"
print(os.path.exists(rd), sorted(os.listdir(rd)))
print(open("step_01.py").read() if os.path.exists("step_01.py") else "no step_01")
