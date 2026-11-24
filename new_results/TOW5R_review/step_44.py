import subprocess, re
t5 = '<project>/research/08_theory/TOW5_tower_depth10'
# 1) diff normalized_summand + total between prove_S34.py and prove_S5.py
def extract(fname, func):
    s = open(t5+'/'+fname).read()
    i = s.find('def '+func)
    j = s.find('\ndef ', i+1)
    return s[i:j]
for func in ['normalized_summand','total','ff','tf','tf1']:
    a = extract('prove_S34.py', func); b = extract('prove_S5.py', func)
    print(func, 'IDENTICAL' if a==b else 'DIFFERENT')
# 2) where does "67/67" come from?
r = subprocess.run(['grep','-l','67/67', '-r', t5], capture_output=True, text=True)
print('67/67 in:', r.stdout)
r = subprocess.run(['grep','-rn','d in {12, 15, 18, 22, 26, 30}', t5, '--include=*.py'], capture_output=True, text=True)
print(r.stdout[:800])
