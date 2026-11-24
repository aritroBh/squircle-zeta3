# Determinants of the per-weight matrices (words ending in 1 against monomials in the Lyndon generators)
# built by g9_hex5_s2_atoms.py (its own code, executed up to the inversion step; nothing is written to disk).
import os, tempfile
os.environ.setdefault('G9WORK', tempfile.gettempdir())
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'g9_hex5_s2_atoms.py')).read()
src = src.split('@lru_cache(None)\ndef reduce_word')[0]
ns = {}
exec(compile(src, 'g9_hex5_s2_atoms_head', 'exec'), ns)
import sympy as sp
for w in range(1, 6):
    words, mons, Minv = ns['BASIS_MAT'][w]
    print(f"weight {w}: {len(words)} words, {len(mons)} monomials, det(M) = {1/Minv.det()}, M invertible over Q: True")
