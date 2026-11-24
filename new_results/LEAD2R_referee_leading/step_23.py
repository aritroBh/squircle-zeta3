import os
print(os.listdir('../TOW2_tower_all_depths'))
try:
    t = open('../TOW2_tower_all_depths/result.md').read()
    print(t[:2500])
except Exception as e:
    print(e)
