import pickle
H4='<project>/research/08_theory/A3h4_frame_ds'
for n in (2,3,4):
    d = pickle.load(open(f'{H4}/validated_w{n}.pickle','rb'))
    print('w%d:'%n, type(d), list(d.keys()) if isinstance(d,dict) else len(d))
    rows = d['rows']
    print('  rows:', len(rows), 'first row type:', type(rows[0]))
    r0 = rows[0]
    if isinstance(r0, dict):
        items = list(r0.items())[:3]
        print('  sample items:', items)
    elif isinstance(r0, (list,tuple)):
        print('  sample:', r0[:3])
