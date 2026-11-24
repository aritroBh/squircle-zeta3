base = '<project>/research/08_theory'
c = base + '/ZQ5c_zagier_w5_capped'
print(open(c+'/summary.md').read()[:3000])
print('===== progress.log tail')
print('\n'.join(open(c+'/progress.log').read().splitlines()[-15:]))
