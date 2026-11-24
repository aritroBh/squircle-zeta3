import glob, os
hits = sorted(glob.glob('../../08_theory/D17W_draft17_writer/mirror_backup/*.tex'))
for h in hits: print(h)
hits2 = sorted(glob.glob('../../08_theory/D17W_draft17_writer/mirror_backup/**/*.tex', recursive=True))
print(len(hits2))
for h in hits2[:60]: print(h)
