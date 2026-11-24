src = open('../../08_theory/D17W_draft17_writer/mirror_backup/main.tex').read()
print(len(src))
# search for P'''(1) and non-local, N_3
import re
for pat in ["P'''", "non-local", "N_3", "N3", "twelve edges", "Future", "future"]:
    idxs = [m.start() for m in re.finditer(re.escape(pat), src)]
    print(pat, len(idxs), idxs[:10])
