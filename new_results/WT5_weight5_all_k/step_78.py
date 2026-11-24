import json
W='<project>/research/08_theory/WT5_weight5_all_k/'
out = {}
cb = json.load(open(W+'cobracket_e.json'))
for nm, cf in [('e1','e1_closed_form.json'), ('e2','e2_closed_form.json'), ('e3','e3_closed_form_rel0.json')]:
    d = json.load(open(W+cf))
    out[nm.upper()] = dict(
        pairing=json.load(open(W+'pairings.json'))[nm.upper()]['pairings'],
        depth_two=cb[nm.upper()]['correction'],
        cobracket_target=cb[nm.upper()]['target'],
        classical_R=dict(names=d['names'], ca=d['ca'], cb=d['cb'], aL=d['aL'], bL=d['bL']),
        stats=d['stats'], relation_resids=d['relation_resids'],
        n_identical_relations=d['n_identical_relations'])
# E3 gamma
out['E3']['gamma_lll'] = json.load(open(W+'e3_gamma_lll.json'))
out['E3']['closed_form_text'] = 'E3 = lam^5/48 + 5*pi^2*lam^3/144 + pi^4*lam/64 + 5*zeta(3)*lam^2/16 + 5*zeta(3)*pi^2/96 + 15*zeta(5)/32'
# W5pp reference
w5cf = json.load(open(W+'w5pp_closed_form.json'))
n5 = json.load(open(W+'cobracket5pp.json'))
out['W5pp'] = dict(depth_two=n5['correction'], cobracket_target=n5['target'],
                   classical_R=dict(names=w5cf['names'], ca=w5cf['ca'], cb=w5cf['cb'], aL=w5cf['aL'], bL=w5cf['bL']),
                   stats=w5cf['stats'], relation_resids=w5cf['relation_resids'])
json.dump(out, open(W+'closed_forms.json','w'), indent=1)
print('closed_forms.json written;', {k: len(v['classical_R']['names']) for k,v in out.items()})
