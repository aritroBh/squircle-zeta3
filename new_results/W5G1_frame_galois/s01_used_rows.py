"""s01: extract the used unfolded rows (family, index, label) and their per-sector y from W5C2 certificates.
Writes used_rows.pkl: {tag: [(ufam, uidx, lam(str), y_rat(str), y_sq(str))]}"""
import json, pickle, collections
W5C2='<project>/research/08_theory/W5C2_exact_certificate'
W5C1=W5C2+'/../W5C1_certificate_scoping'
out={}
for tag in ('pp','pm','mp','mm'):
    cert=json.load(open(f'{W5C2}/certificate_{tag}.json'))
    L=[]
    for ent in cert['rows']:
        s=ent['src']
        if isinstance(s,str): s=eval(s)
        L.append((s['ufam'],int(s['uidx']),str(s['lam']),ent['y_rat'],ent['y_sq']))
    out[tag]=L
    print(tag,len(L),collections.Counter(x[0] for x in L))
pickle.dump(out,open('used_rows.pkl','wb'))
lift1=pickle.load(open(W5C1+'/rels_lift1.pkl','rb'))
for fam in ('lift1',):
    c=collections.Counter()
    for tag,L in out.items():
        for f,i,*_ in L:
            if f==fam: c[(tag,lift1[i][0][1])]+=1
    print(fam,'by base index',sorted(c.items()))
