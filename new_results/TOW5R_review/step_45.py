t5 = '<project>/research/08_theory/TOW5_tower_depth10'
# show the tf1 and total diffs
def extract(fname, func):
    s = open(t5+'/'+fname).read()
    i = s.find('def '+func)
    j = s.find('\ndef ', i+1)
    return s[i:j]
print('S34 tf1:', repr(extract('prove_S34.py','tf1')))
print('S5  tf1:', repr(extract('prove_S5.py','tf1')))
print('S34 total:', repr(extract('prove_S34.py','total')))
print('S5  total:', repr(extract('prove_S5.py','total')))
