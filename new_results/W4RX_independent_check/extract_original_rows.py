"""W4RX: RERUN of the ORIGINAL code (the only W4RX script that imports modules of the notes).

Purpose: produce, with the original pipeline exactly as verify_D.py uses it, the 5068 rows named
by D_certificate.json, the full list of original weight-4 standard rows, and the original target
D, so that the independent scripts can COMPARE against them. Nothing produced here is used to
build anything in the independent scripts. All logging of the original modules is redirected to
stdout here (the original log files are not touched).
"""
import sys, json, os
sys.dont_write_bytecode = True
ROOT = '<project>/research/08_theory'
sys.path.insert(0, ROOT + '/A3h5_exact_certificate')
import verify_D as V          # imports frame_ds, cmzv, relations through its own sys.path lines
def mylog(s): print('[original code] ' + s, flush=True)
V.log = mylog
import relations, cmzv, frame_ds
relations.log = mylog; cmzv.log = mylog; frame_ds.log = mylog

cert = json.load(open(ROOT + '/A3h5_exact_certificate/D_certificate.json'))
std, bases = V.rebuild_standard()
out_rows = []
for label in cert['labels']:
    if label[0] == 'T0':
        lab = label[1]
        if lab[0] == 'std':
            row = std[lab[1]]
        else:
            n = {'lift1': 1, 'lift2old': 2, 'lift3old': 3}[lab[0]]
            row = V.primitive(V.mul(dict(bases[n][lab[1]]), {tuple(lab[2]): 1}))
    else:
        _, n, i, w = label
        data = V.load(V.HERE / f'validated_w{n}.pickle')
        g, fs = data['provenance'][i]
        p0, p1, pinf = g
        ts = {V.cross(p, p0, p1, pinf): p for p in V.POINTS}
        raw = V.primitive(V.convert(g, ts, V.identity(fs)))
        assert raw == data['rows'][i]
        row = V.primitive(V.mul(dict(raw), {tuple(w): 1}))
    out_rows.append([[list(w), str(c)] for w, c in row])
targets = V.rebuild_target()
json.dump({'labels': cert['labels'],
           'rows': out_rows,
           'std_all': [[[list(w), str(c)] for w, c in row] for row in std],
           'bases': {str(n): [[[list(w), str(c)] for w, c in row] for row in bases[n]] for n in bases},
           'D': [[[list(w), str(c)] for w, c in sorted(t.items())] for t in targets]},
          open('data/original_rows.json', 'w'))
print('wrote data/original_rows.json', len(out_rows), 'rows;', len(std), 'std rows')
