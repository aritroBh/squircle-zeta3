from common import *
for W in (5,6,7,8):
    data,B = load_rows_meta(W)
    bad=0; n=0
    for fam,tag,_ in data:
        r = my_relation(fam,tag); n+=1
        if not all(conv(w) and len(w)==W for w in r): bad+=1
    print(f"W={W}: {n} relations regenerated from labels; relations containing a divergent/wrong-weight word: {bad}")
