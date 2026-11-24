import os
m5b='<project>/research/08_theory/M5b_w6_closed_form'
src=open(os.path.join(m5b,'w6_values.py')).read()
old = src[src.index('def main(digits):'):src.index("if __name__ == '__main__':")]
new = '''def main(digits):
    ctx.prec = math.ceil((digits+140)*math.log2(10))
    mp.mp.dps = digits+120
    cache = HERE/f'w6_leaves_{digits}.pkl'
    if cache.exists():
        leaves = {eval(k): mp.mpf(v) for k, v in pickle.loads(cache.read_bytes()).items()}
        log(f'loaded leaf cache ({len(leaves)} leaves)')
    else:
        t0 = time.monotonic()
        leaves = {}
        for rec in list(get_leaves(digits)):
            v = leaf_value(rec, digits+40)
            leaves[rec] = mp.mpf(v.mid().str(digits+110, radius=False))
            log(f'  leaf {rec} done ({time.monotonic()-t0:.1f}s)')
        cache.write_bytes(pickle.dumps({repr(k): mp.nstr(v, digits+110) for k, v in leaves.items()}))
        log(f'leaf values computed ({time.monotonic()-t0:.1f}s)')
    # atom values for the full basis
    basis = atoms_w6(BASE_LI6)
    acache = HERE/f'w6_atomvals_{digits}.pkl'
    if acache.exists():
        vals = [mp.mpf(v) for v in pickle.loads(acache.read_bytes())]
        log(f'loaded atom value cache ({len(vals)})')
    else:
        t0 = time.monotonic()
        vals = []
        for name, _, recipe in basis:
            if recipe[0] == 'product':
                v = mp.mpf(1)
                for leaf in recipe[1]:
                    v *= leaves[leaf]
                vals.append(v)
            elif recipe[0] in ('li', 'mpl'):
                vals.append(leaves[recipe])
            else:
                raise ValueError(recipe)
        acache.write_bytes(pickle.dumps([mp.nstr(v, digits+110) for v in vals]))
        log(f'atom values computed ({time.monotonic()-t0:.1f}s)')
    return leaves, vals

'''
src = src.replace(old, new)
# rewrite the control block (uses mpf arithmetic now)
old2 = src[src.index("if __name__ == '__main__':"):]
new2 = '''if __name__ == '__main__':
    digits = int(sys.argv[1]) if len(sys.argv) > 1 else 1650
    leaves, vals = main(digits)
    d1 = leaves[('li', 6, R/2)] + leaves[('li', 6, -R/2)] - leaves[('li', 6, Q(F(1,2)))]/32
    print('Li6 duplication control (1/sqrt2 pair vs 1/2):', mp.nstr(abs(d1), 5))
    d2 = leaves[('li', 6, B)] + leaves[('li', 6, -B)] - leaves[('li', 6, B*B)]/32
    print('Li6 duplication control (b pair vs b^2):', mp.nstr(abs(d2), 5))
    import w5_search as WS5
    ctx.prec = math.ceil((digits+140)*math.log2(10))
    ref = WS5.mpl(3, 2, B, ONE, digits)
    dd = leaves[('mpl', 3, 2, B, ONE)] - mp.mpf(ref.mid().str(digits+110, radius=False))
    print('Li32(b,1) mpl_general vs M4c mpl:', mp.nstr(abs(dd), 5))
    n6v = json.loads((HERE.parent/'M5a_next_coefficient/n6_values.json').read_text())
    mp.mp.dps = 120
    for key, rec in [('Li33_b1', ('mpl', 3, 3, B, ONE)), ('Li33_s1', ('mpl', 3, 3, R/2, ONE)),
                     ('Li33_b2', ('mpl', 3, 3, B, Q(2))), ('Li42_bs', ('mpl', 4, 2, B, R/2))]:
        print(key, 'vs M5a n6_values:', mp.nstr(abs(leaves[rec] - mp.mpf(n6v[key])), 5))
    mp.mp.dps = 80
    zz = Q(F(1,2))-R/4
    zf = mp.mpf(1)/2 - mp.sqrt(2)/4
    direct = mp.polylog(6, zf)
    print('Li6((2-sqrt2)/4) acb vs mpmath polylog:', mp.nstr(abs(leaves[('li', 6, zz)]-direct), 5))
'''
src = src.replace(old2, new2)
open(os.path.join(m5b,'w6_values.py'),'w').write(src)
import py_compile; py_compile.compile(os.path.join(m5b,'w6_values.py'), doraise=True)
print('ok')
import subprocess
PY='<project>/.venv/bin/python'
r=subprocess.run(['nice','-n','10',PY,'-B','w6_values.py','1650'],capture_output=True,text=True,cwd=m5b,timeout=1700)
print(r.stdout[-2500:]); print(r.stderr[-900:])
