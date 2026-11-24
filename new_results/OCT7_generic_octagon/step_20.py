import os
base="research/08_theory"
# Look at D15b_polygon_law to see how term values are computed for polygons
p=os.path.join(base,"D15b_polygon_law")
print(sorted(os.listdir(p)))
print("="*70)
print(open(os.path.join(p,"result.md")).read()[:5000])
