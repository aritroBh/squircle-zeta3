import os, json
d17='research/08_theory/D17_zeta5_polygons'
tr=json.load(open(os.path.join(d17,'transcript.json')))
# print the last messages of the transcript
if isinstance(tr, list):
    for m in tr[-3:]:
        s=str(m)
        print(s[:4000])
        print('-'*80)
else:
    print(str(tr)[:4000])
