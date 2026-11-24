
import os
p='research/10_frontier/E1_edge_law/verify'
print(os.listdir(p))
src=open(os.path.join(p,'vjet.py')).read()
print(src[:3000])
