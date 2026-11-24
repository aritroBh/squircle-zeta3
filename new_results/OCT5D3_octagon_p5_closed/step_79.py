import os, sys, pickle, time
os.chdir("<project>/research/08_theory/OCT5D3_octagon_p5_closed")
sys.path.insert(0, ".")
import mpmath as mp
mp.mp.dps = 1650
t0 = time.time()
import importlib, oct5d3_atoms
importlib.reload(oct5d3_atoms)
A = oct5d3_atoms.base_atoms()
print("recomputed base_atoms in %.1fs, %d atoms" % (time.time()-t0, len(A)))
print([k for k in A if 'Li4' in k])
pickle.dump({k: mp.nstr(v, 1630) for k, v in A.items()}, open("atoms_1650.pkl","wb"))
print("saved")
