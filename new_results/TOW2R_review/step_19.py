import json, os, re
tow2 = '<project>/research/08_theory/TOW2_tower_all_depths'
tr = json.load(open(os.path.join(tow2, 'transcript.json')))
txt = json.dumps(tr)
for pat in [r'V\(3\) non-adjacent remainder: [^\\]*', r'V\(4\) non-adjacent remainder: [^\\]*',
            r'E4-slope \(4,1\) non-adjacent remainder: [^\\]*', r'adjacent \(\d,\d\) remainder: \d+']:
    ms = re.findall(pat, txt)
    print(pat, '->', ms[:8])
