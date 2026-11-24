"""Which rows do the 74 certificates reference (from the core files' frame and t sections)? -> used_rows.json"""
from common import *
used = {'framew5': set(), 'frame234': set(), 'frame234_base': set(), 'A3h4_frame_base': set(), 'targets_row': set(), 'zero_coeffs': 0}
for name in TARGETS:
    for o in read_jsonl_xz(f'{STORE}/{name}.core.jsonl.xz'):
        if o['section'] == 'meta': used['targets_row'].add(o['i'])
        if o['section'] == 'frame':
            for t in o['terms']:
                if rat(t[-1])[0] == 0: used['zero_coeffs'] += 1; continue
                if t[0] == 'framew5': used['framew5'].add(t[1])
                else: used['frame234'].add(t[1]); used['frame234_base'].add((t[2], t[3]))
        if o['section'] == 't':
            for t in o['terms']:
                if t[0] == 'A3h4_frame' and rat(t[-1])[0] != 0: used['A3h4_frame_base'].add((t[1], t[2]))
out = {k: (sorted(v) if isinstance(v, set) else v) for k, v in used.items()}
json.dump(out, open(HERE + '/used_rows.json', 'w'))
log('used rows: ' + json.dumps({k: (len(v) if isinstance(v, list) else v) for k, v in out.items()}) + ' frame234 bases ' + str(out['frame234_base']) + ' A3h4 bases ' + str(out['A3h4_frame_base']), 'used_rows.log')
