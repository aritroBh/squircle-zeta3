import sys, time, pickle, os
work = '<project>/research/08_theory/TOW5_tower_depth10'
sys.path.insert(0, work)
import trunc_engine as te
t0=time.time()
S2 = te.stirling_table(45)
RHO = {}
for d in range(7, 41):
    rho_d = te.rho_level(d, 10, S2)
    RHO[d] = {key: [tuple([c for c in (v.p0,v.p1,v.p2,v.p3)]) for v in vals] for key, vals in rho_d.items()}
pickle.dump(RHO, open(work+'/data_rho_D10_d7_40.pkl','wb'))
print(f"rho data d=7..40 D=10 generated ({round(time.time()-t0,1)}s)")
# count cells for rho(9), rho(10)
n9 = sum(len(RHO[d]) for d in range(7,41) if 2*d-4 >= 9)
n10 = sum(len(RHO[d]) for d in range(7,41) if 2*d-4 >= 10)
print("cells: rho(9):", n9, " rho(10):", n10)
with open(work+'/progress.log','a') as f:
    f.write(f"{time.strftime('%H:%M:%S')} stepB1: rho data d=7..40 D=10 generated; cells rho9={n9} rho10={n10}\n")
