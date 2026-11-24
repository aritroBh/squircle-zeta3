base = '<project>'
print(open('step_01.py').read() if os.path.exists('step_01.py') else 'no step_01')
print('='*100)
print(open(base + '/research/08_theory/G2_M4alpha/result.md').read()[:6000])
