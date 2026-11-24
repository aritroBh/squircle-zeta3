import os

base = "<project>/research"
cor4 = os.path.join(base, "08_theory/COR4_M3alpha_classical")
print("COR4 dir:", os.listdir(cor4))
print()
cwd = os.getcwd()
print("CWD:", cwd)
print("CWD contents:", os.listdir(cwd) if os.path.exists(cwd) else "MISSING")
