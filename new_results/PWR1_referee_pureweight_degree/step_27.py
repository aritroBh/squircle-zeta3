import os
print([f for f in os.listdir('../../10_frontier/H1_octahedron_end') if f.endswith(('.py','.out','.txt','.json'))])
