import os, glob
print(os.getcwd())
for f in sorted(glob.glob('*')) + sorted(glob.glob('../D11b_stokes_numeric/*')):
    print(f)
