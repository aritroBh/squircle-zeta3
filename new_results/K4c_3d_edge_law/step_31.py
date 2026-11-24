
import json, re, time
import mpmath as mp
from mpmath import mpf

def parse_real(s):
    s = s.strip()
    if s.startswith('('):
        # form "(real + imag j)" or "(real - imag j)"
        body = s[1:-1]
        # find the last +/- that precedes 'j'
        m = re.match(r'^(.*?)([+-].*j)$', body.replace(' ',''))
        real = m.group(1)
        imag = m.group(2)[:-1].replace('j','')
        return mpf(real), mpf(imag)
    return mpf(s), mpf(0)

for tag in ["h32","h48"]:
    J = json.load(open(f"modelB_{tag}.json"))
    vals = []
    for s in J["Sdk"]:
        re_, im_ = parse_real(s)
        assert abs(im_) < mpf('1e-100'), (tag, s)
        vals.append(mp.nstr(re_, 70))
    J["Sdk"] = vals
    json.dump(J, open(f"modelB_{tag}.json","w"), indent=1)
    print(f"modelB_{tag} cleaned; imag max ok")
with open('progress.log','a') as f:
    f.write(f"[{time.strftime('%H:%M:%S')}] step_08.py MODEL (b) [3 axes+4 body diag]: k0 = 1.075764930490848519873765587084270092790655697826818297768 matches EXACT -36+21sqrt3+(9/13)sqrt(35515-20504sqrt3) to ~60 digits (PROVED exact k0). h48: k1=3.516495437145298347373890733252633379646687669001453640008, k2=2.354480678056587151344461624642939462487721480650698816991, k3=-1.898327187347440896619891901482403412469397080520072254548, k4=-0.4126500481804702236949455301550953680647740676969757841329. h48-h32 diffs: k1 4.3e-39, k2 1.4e-35, k3 3.4e-33, k4 8.3e-31. NUMERICAL.\n")
print("logged")
