import subprocess, time
cwd='<project>/research/08_theory/DCYA1_corner_decay_proof/'
with open(cwd+'progress.log','a') as f:
    f.write("step: VERIFIED exact representation M(u;a)=T(u)+R(u) to ~50 digits at (u,a)=(0.1,0.3)x(pi/3,2pi/3,3pi/4,pi/2). T=int e^{-(1+u)L}e^{-x} Btilde((u-1)x) dx, R=int[e^{-(1+u)L}-1]e^{-ux}dx, Btilde(z)=(e^z-2cos a)/(1+sqrt(1-2cos a e^z+e^{2z})).\n")
    f.write("step: built accurate M_k tables (custom series+tanh-sinh, high dps): pi/2 k<=100 (dps140, matches K2 M_100), 2pi/3 & pi/3 k<=80 (dps120). FFT method has noise floor at large k, abandoned.\n")
    f.write("step: main term identified = K2 Hankel at z=+-i alpha, m=1 weight, reduced to 1D integrals I_+^{(j)} (Tricomi-U form). leading+g_1 captures kappa_1 (rest ~ N^{-2}e^{-theta}). Verified numerically.\n")
print(open(cwd+'progress.log').read())
