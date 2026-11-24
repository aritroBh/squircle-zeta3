import json, os
d = json.load(open('../E5d_nonlocal_pslq/b0_p5_convB_350.json'))
e1, e2 = d['e1']['P5'], d['e2']['P5']
print("len e1:", len(e1), "len e2:", len(e2))
# count digits after decimal
i1 = e1.index('.'); 
d1 = len(e1) - i1 - 1
i2 = e2.index('.')
d2 = len(e2) - i2 - 1
print("digits after decimal: e1 =", d1, ", e2 =", d2)
# common prefix
n = min(len(e1), len(e2))
k = 0
while k < n and e1[k] == e2[k]:
    k += 1
print("common prefix length (chars):", k)
print("agreement digits after decimal:", k - i1 - 1)
print("first diff: e1[%d]=%r e2[%d]=%r" % (k, e1[k-3:k+3], k, e2[k-3:k+3]))
print("conv estimates:", d['e1']['conv'], d['e2']['conv'])
# existing step files
for f in ['step_01.py','step_02.py']:
    if os.path.exists(f):
        print("==== ", f)
        print(open(f).read()[:2000])
