import os
p = '../TOWN1_tower_large_n'
print(os.listdir(p))
t = open(p+'/result.md').read()
print(t[:5500])
