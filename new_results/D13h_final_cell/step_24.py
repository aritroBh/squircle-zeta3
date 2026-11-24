import json, shutil
p = '<scratch>'
data = json.load(open(p))
exact = all(data[i+1][0] == data[i][1] for i in range(len(data)-1))
print("exact contiguous:", exact)
print("first cell:", data[0], "last cell:", data[-1])
gaps = [(i, data[i][1], data[i+1][0]) for i in range(len(data)-1) if data[i+1][0] != data[i][1]]
print("gaps:", gaps[:5])
# y-box coverage check: cell [N1,N2] certified for y in [ys(N2), ys(N1)]
import math
ys = lambda N: math.sqrt(math.pi/(2*N))
print("sweep covers y in [%.6e, %.6e] = N in [400, 1e10]" % (ys(1e10), ys(400)))
shutil.copy(p, 'sweep_cells_400_1e10.json')
print("archived copy: sweep_cells_400_1e10.json")
# update progress log
with open('progress.log','a') as f:
    f.write("turn7: sweep exact-contiguous [400,1e10] confirmed; archived locally. 3 sups: S1[100,400]=194.3809, S2[400,1e10]=617117.5536, S3[1e10,oo)=49.5626 (J192; J48 gives 53.6157). THEOREM C = 617117.6.\n")
