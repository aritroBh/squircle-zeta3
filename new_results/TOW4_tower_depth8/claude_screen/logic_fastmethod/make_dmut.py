slow = open('prove_adj56.py').read()
fast = open('prove_adj56_fast_T3.py').read()
patches = {
 'fold_sign':      ('terms.append(-p * sp.expand(g)', 'terms.append(p * sp.expand(g)'),
 'fold_phase':     ('i**(-6-delta)', 'i**(-5-delta)'),
 'fold_mu':        ('mu = 4*p-4-t', 'mu = 4*p-5-t'),
 'rA0':            ('rA0 = i*(4*p-3)/(p+1)', 'rA0 = i*(4*p-3)/(p+2)'),
 'rB0':            ('rB0 = (4*p-1)*(4*p-3)/((p+1)*p)', 'rB0 = (4*p+1)*(4*p-3)/((p+1)*p)'),
 'lhs_s':          ('eR(delta, s+1, d, p)', 'eR(delta, s, d, p)'),
 'Ediag_s_shift':  ('Ediag(dp, s)/(2*d-dp-1))', 'Ediag(dp, s+1)/(2*d-dp-1))'),
 'Ediag_altdrop':  ('tot += (-1)**t * ehat(delta, s+t, conj)', 'tot += ehat(delta, s+t, conj)'),
 'conj_tau1_drop': ('(-1)**n * T3.tau2(n) if conj else T3.tau1(n)', 'T3.tau2(n) if conj else T3.tau1(n)'),
 'binom_shift':    ('binom = ff(mu-2, delta-t)/sp.factorial(delta-t)', 'binom = ff(mu-1, delta-t)/sp.factorial(delta-t)'),
 'family1_tau':    ('terms.append((p+1)*tau1(delta-dp)*eR(dp, s, 2*p, p+1)*rA0/(2*d-dp-1))', 'terms.append((p+1)*tau2(delta-dp)*eR(dp, s, 2*p, p+1)*rA0/(2*d-dp-1))'),
}
for k,(o,n) in patches.items():
    for tag,src in (('slow',slow),('fast',fast)):
        assert src.count(o)==1, (k,tag,src.count(o))
        open(f'dmut/{tag}_{k}.py','w').write(src.replace(o,n))
print(len(patches)*2,'files')
