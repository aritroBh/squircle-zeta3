import os, subprocess
files = sorted([f for f in os.listdir('.') if f.startswith('step_') and f.endswith('.py')])
print(len(files), 'step files')
# identify key ones
for f in files:
    head = open(f).read(300)
    tag = ''
    if 'Mobius maps' in head or 'mobius_from_triple' in head: tag = 'MAPS'
    if 'cross-ratio' in head or 'permutations(PTS,4)' in head: tag += ' XRAT'
    if 'cert_w23' in head: tag += ' CERT23'
    if 'w4_y.pkl' in head: tag += ' CERT4'
    if 'verify_hexc2' in head: tag += ' VERIF'
    if 'result.md' in head or 'result = ' in head[:50]: tag += ' RESULT'
    if tag: print(f, head.splitlines()[0][:60], '=>', tag)
