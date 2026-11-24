import os, json, subprocess
os.chdir('<project>/research/08_theory/W5C2_exact_certificate')
with open('progress.log','a') as f:
    f.write('2026-10-29 06:32:14 step_05: EXACT UNFOLDED IDENTITY HOLDS for D_rat and D_sq on ALL 46656 words (Fractions, residual 0), from certificates + source rows.\n')
    f.write('2026-10-29 06:35:36 step_06: Galois transforms g.r of used rows: 36016/40392 literally in row files (std all present); missing 4376 (lift1/lift2old/lift4old/framew5 transforms) - proof does NOT need set-closure ((g.r).V = g(r.V) = 0 via Galois conjugation of values).\n')
for tag in ['pp','pm','mp','mm']:
    sz = os.path.getsize(f'certificate_{tag}.json')
    cert = json.load(open(f'certificate_{tag}.json'))
    r0 = cert['rows'][0]; r1 = cert['rows'][1]
    print(tag, sz//10**6, 'MB; rows', len(cert['rows']), 'rank', cert['rank'], 'fam', cert['family_counts_used'])
    print('  row0:', r0['fam'], r0['famidx'], 'src', r0['src'], 'y_rat', r0['y_rat'][:30], 'tuple[:3]', r0['tuple'][:3])
print(subprocess.run(['bash','-c','ls -la certificate_*.json work/*.pkl | head -30'], capture_output=True, text=True).stdout)