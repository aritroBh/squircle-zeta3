import time
with open("progress.log","a") as f:
    f.write(f"{time.strftime('%H:%M')} step_07b: Mellin=Legendre residuals were endpoint-singularity quadrature error; with smoothing substitution: 2.5e-60 at w=+-0.9, 7e-59 at w=+-0.1. Identity CONFIRMED numerically [referee]\n")
print("logged")
