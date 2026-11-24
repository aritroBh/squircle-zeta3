import os
def find_files_with(d, needle, ext='.py', maxn=10):
    out=[]
    for root,dirs,files in os.walk(d):
        for f in files:
            if f.endswith(ext):
                p=os.path.join(root,f)
                try:
                    if needle in open(p, errors='ignore').read():
                        out.append(p)
                except: pass
    return out[:maxn]

# In K1_nonlocal_law, find where kappa and A0,A1,B are set from an edge
hits = find_files_with('research/08_theory/K1_nonlocal_law','kappa')
print("K1_nonlocal_law files with 'kappa':", [os.path.basename(h) for h in hits])
# Look for the parametrization comment "A1+B=1" or "B(1+"
for d in ['K1_nonlocal_law','K1c_engine_predict','K1d_engine_fixed']:
    dd='research/08_theory/%s/'%d
    hits2 = find_files_with(dd,'A1+B')
    print(d, "files with 'A1+B':", [os.path.basename(h) for h in hits2])
