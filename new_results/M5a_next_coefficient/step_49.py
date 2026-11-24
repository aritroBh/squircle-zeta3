import subprocess, os
work = '<project>/research/08_theory/M5a_next_coefficient'
# Identify auto-saved step files by content
for n in range(18, 49):
    f = f'{work}/step_{n:02d}.py'
    if os.path.exists(f):
        src = open(f).read()
        first = [l for l in src.splitlines() if l.strip() and not l.startswith('import')][:1]
        tag = ''
        if 'omega5.json' in src and 'reconstruction' in src: tag = 'OMEGA5 PARSE'
        if 'results_summary' in src: tag = 'SUMMARY'
        if 'witt' in src.lower(): tag = 'LIE DIMS'
        if 'li_mn' in src: tag = 'ATOM SERIES'
        if 'W3exact' in src: tag = 'TAILS CHECK'
        if 'real_g_raw' in src and 'n6_values' in src: tag = 'ATOM INTEGRAL XCHECK'
        print(n, tag, first[:1])