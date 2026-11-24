import os
os.chdir("<project>")
os.makedirs("research/08_theory/SRF2_cube_residue", exist_ok=True)
os.chdir("research/08_theory/SRF2_cube_residue")
print("cwd:", os.getcwd())
print(os.listdir("."))
import sympy, mpmath, numpy, scipy
print("sympy", sympy.__version__, "| mpmath", mpmath.__version__, "| numpy", numpy.__version__)
try:
    import flint
    print("flint ok")
except Exception as e:
    print("flint:", e)
