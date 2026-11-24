
import re
def parse_complex(s):
    s = s.strip()
    # find last +/- not preceded by e/E
    idx = None
    for i in range(len(s)-1, 0, -1):
        if s[i] in '+-' and s[i-1] not in 'eE':
            idx = i; break
    re_part, im_part = s[:idx], s[idx:]
    im_part = im_part.replace('j','').replace('i','')
    return re_part, im_part

def load(fname):
    d = {}
    with open(fname) as f:
        for line in f:
            line=line.strip()
            if not line or line.startswith('#'): continue
            parts = line.split('\t')
            if len(parts)<2: parts = line.split()
            k = int(parts[0])
            r,i = parse_complex(parts[1])
            d[k] = (r,i)
    return d

base='research/08_theory/K2_corner_decay/out/'
files = {'pi_3':'mk_pi_3.tsv','pi_2':'mk_pi_2.tsv','2pi_3':'mk_2pi_3.tsv','3pi_4':'mk_3pi_4.tsv','5pi_6':'mk_5pi_6.tsv'}
data = {a: load(base+f) for a,f in files.items()}
for a,d in data.items():
    ks = sorted(d)
    print(a, 'k range', ks[0], ks[-1], 'n=',len(ks))
    for k in [0, 100, 300]:
        print('  k=%d re=%s im=%s'%(k, d[k][0][:45], d[k][1][:45]))
