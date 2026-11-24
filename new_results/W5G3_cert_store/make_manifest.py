import os, hashlib
rows = []
for d in ('store',):
    for n in sorted(os.listdir(d)):
        p = os.path.join(d, n); h = hashlib.sha256()
        with open(p, 'rb') as f:
            for b in iter(lambda: f.read(1 << 22), b''): h.update(b)
        rows.append((p, os.path.getsize(p), h.hexdigest()))
with open('MANIFEST.txt', 'w') as f:
    f.write('# file  size_bytes  sha256   (xz files; core = frame a_j + W + t_k, lit = literal row coefficients)\n')
    for r in rows: f.write(f'{r[0]}  {r[1]}  {r[2]}\n')
core = sum(r[1] for r in rows if '.core.' in r[0]); lit = sum(r[1] for r in rows if '.lit.' in r[0])
print(f'{len(rows)} files; core {core/1e9:.3f} GB; lit {lit/1e9:.3f} GB; total {(core+lit)/1e9:.3f} GB')
