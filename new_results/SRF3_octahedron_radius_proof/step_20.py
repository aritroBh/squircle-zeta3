
import os
base = "<project>"
# tail of the K80 run output for timing
p = os.path.join(base,"research/10_frontier/H1_octahedron_end/hi_n3_K80.out")
src = open(p).read()
print(src[-800:])
print("SIZE:", os.path.getsize(p))
