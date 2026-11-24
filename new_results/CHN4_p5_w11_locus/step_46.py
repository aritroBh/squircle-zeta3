import numpy as np
# res from previous cell: (a1,a2,a4,a5,q)
arr=np.array([(r[0],r[1],r[2],r[3],r[4]) for r in res])
# Look at log-log map of sign
pos=arr[arr[:,4]>0]; neg=arr[arr[:,4]<0]
print("positive Q11 points: a1 range", pos[:,0].min(), pos[:,0].max(), " a2 range", pos[:,1].min(), pos[:,1].max())
print("negative Q11 points: a1 range", neg[:,0].min(), neg[:,0].max(), " a2 range", neg[:,1].min(), neg[:,1].max())
# sample some positive points with moderate values
mod=pos[np.abs(pos[:,4])<1e3]
print("\nsample moderate positive-Q11 V10+ points (a1,a2,a4,a5,Q11):")
for r in mod[:8]:
    print(f"  a1={r[0]:.4g} a2={r[1]:.4g} a4={r[2]:.4g} a5={r[3]:.4g} Q11={r[4]:.6g}")
print("\nsmallest-|Q11| points (sign-change candidates):")
idx=np.argsort(np.abs(arr[:,4]))[:8]
for r in arr[idx]:
    print(f"  a1={r[0]:.4g} a2={r[1]:.4g} a4={r[2]:.6g} a5={r[3]:.6g} Q11={r[4]:.6g}")
