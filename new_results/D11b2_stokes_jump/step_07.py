import json
t = json.load(open('../D11b_stokes_numeric/transcript.json'))
print(type(t), len(t))
# print the last few entries' text briefly
def show(e, n=1500):
    s = str(e)
    print(s[:n]); print('---')
for e in t[-6:]:
    show(e)
