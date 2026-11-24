import os
print("CWD:", os.getcwd())
base='research/10_frontier/E1_edge_law'
print(open(os.path.join(base,'s01_out.txt')).read())
print("="*70)
# check s13_general_m.py briefly (regular 2m-gon edge data)
src=open(os.path.join(base,'s13_general_m.py')).read()
print(src[:3000])
