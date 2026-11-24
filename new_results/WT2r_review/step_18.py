import shutil, os, sys
src = '<project>/research/08_theory/WT2b_stabilization'
work = '<project>/research/08_theory/WT2r_review'
for f in ['fastred.py','genintegrand.py','buildintegrand.py','record_W2b.pkl']:
    shutil.copy(f'{src}/{f}', f'{work}/{f}')
print("copied:", os.listdir(work))
# sanity: load and check record structure
os.chdir(work)
sys.path.insert(0, work)
import pickle
record = pickle.load(open(f'{work}/record_W2b.pkl','rb'))
ks = sorted(set(k for k,_ in record))
print("levels in record:", ks)
print("n records:", len(record))
print("sample keys at k=2:", [kk for kk in record if kk[0]==2])
print("fields:", list(record[(2,(2,0))].keys()))
