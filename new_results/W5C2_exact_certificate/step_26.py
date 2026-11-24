import pickle
for fn in ['rels_std.pkl','rels_frame234.pkl']:
    rows = pickle.load(open(f'../W5C1_certificate_scoping/{fn}','rb'))
    print(fn, len(rows), type(rows[0]))
    lab, row = rows[0]
    print('  label:', lab if len(str(lab))<120 else str(lab)[:120])
    print('  row head:', row[:3], 'types:', type(row[0][0]), type(row[0][1]))
    mx = max(abs(a) for _,r in rows[:50] for _,a in r)
    print('  sample max|cof|:', mx)
vw = pickle.load(open('../W5C1_certificate_scoping/validated_w5.pkl','rb'))
print('validated_w5:', type(vw), list(vw.keys()) if isinstance(vw, dict) else len(vw))