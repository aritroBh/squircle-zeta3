import os, shutil
os.chdir("<project>/research/08_theory/AV3_fractional_fourier")
nested = "research"
if os.path.isdir(nested):
    print(open(os.path.join(nested,"research/08_theory/AV3_fractional_fourier/progress.log")).read() if os.path.exists(os.path.join(nested,"research/08_theory/AV3_fractional_fourier/progress.log")) else "no log")
    shutil.rmtree(nested)
print(os.getcwd(), os.listdir('.'))
for fn in ["step_01.py","step_02.py"]:
    print("="*20, fn)
    print(open(fn).read()[:3000])
