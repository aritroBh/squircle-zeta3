import pickle, os
fd1='<project>/research/08_theory/FD1_angle_factor_in_d'
for f in ['step_01.py','step_02.py']:
    print(f'=== {f} ==='); print(open(fd1+'/'+f).read()[:600])
a2='<project>/research/08_theory/ALPHA2_angle_leading_law'
Fd=pickle.load(open(a2+'/Fd_data.pkl','rb'))
print("Fd_data:", {d: {str(c): v for c,v in sorted(Fd[d].items())} for d in sorted(Fd)})
fits=pickle.load(open(a2+'/Fd_fits.pkl','rb'))
print("Fd_fits keys:", list(fits.keys()) if hasattr(fits,'keys') else fits)
