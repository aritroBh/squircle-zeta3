import os
os.makedirs("research/08_theory/AV3_fractional_fourier", exist_ok=True)
os.chdir("research/08_theory/AV3_fractional_fourier")
print(os.getcwd())
with open("progress.log","a") as f:
    f.write("start: AV3 fractional Fourier path |x|^s, sgn(x)|x|^s\n")
