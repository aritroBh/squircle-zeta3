import sys, pickle; fd1='<project>/research/08_theory/FD1_angle_factor_in_d'
sys.path.insert(0,fd1)
K=pickle.load(open(fd1+'/kern_sq_11_20.pkl','rb'))['kernels']
d=10
keys = [(i, j) for i in range(d, (d-1)//2, -1) for j in [d-i]]
for key in keys:
    print(key, {k: tuple(sorted(K[k][key]['real'].keys())) for k in range(11,21)})
