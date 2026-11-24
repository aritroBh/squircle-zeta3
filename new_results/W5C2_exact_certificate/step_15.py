import pickle, time
t0=time.time()
re_ = pickle.load(open('../W5C1_certificate_scoping/sector/rows_engine_1_1.pkl','rb'))
print("engine pp:", len(re_), type(re_[0]), re_[0][:3])
print("entry types:", type(re_[0][0][0]), type(re_[0][0][1]))
print("max |entry|:", max(abs(a) for r in re_[:2000] for _,a in r), "lens:", min(len(r) for r in re_), max(len(r) for r in re_))
import numpy as np
nnz = np.array([len(r) for r in re_])
print("nnz mean/median:", nnz.mean(), np.median(nnz), "total nnz:", nnz.sum())
rf = pickle.load(open('../W5C1_certificate_scoping/sector/rows_framew5_1_1.pkl','rb'))
print("framew5 pp:", len(rf), type(rf[0]))
r0 = rf[0]
print(list(r0.items())[:3] if isinstance(r0,dict) else r0[:3])
D = pickle.load(open('../W5C1_certificate_scoping/sector/D_rat_1_1.pkl','rb'))
print("D_rat pp:", len(D), type(D), list(D.items())[:2])
print(f"loaded in {time.time()-t0:.1f}s")