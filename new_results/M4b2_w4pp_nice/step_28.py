import json
from fractions import Fraction as F
here = '<project>/research/08_theory/M4b2_w4pp_nice'
cert = json.load(open(here + '/W4pp_nice_certificate.json'))
print('atom count:', len(cert['atoms']), '| labels:', len(cert['labels']))
li4 = [n for n in cert['atom_names'] if n.startswith('Li4')]
print('Li4 atoms:', li4)
cr = [F(a) for a in cert['coefficients'][0]]; cs = [F(a) for a in cert['coefficients'][1]]
atoms_in_cert = [(l[1], cr[i], cs[i]) for i, l in enumerate(cert['labels']) if l[0] == 'atom']
print('atom rows in certificate:', len(atoms_in_cert), '| max |coeff numerator| among Li4:',
      max(max(abs(a.numerator), abs(b.numerator)) for k, a, b in atoms_in_cert if cert['atom_names'][k].startswith('Li4')))
print('last verifier line:', open(here + '/verifier_nice.log').read().strip().splitlines()[-1])
print('progress tail:', open(here + '/progress.log').read().strip().splitlines()[-1])
