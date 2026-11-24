import os
os.chdir('<project>/research/08_theory/LEAD2_leading_proof')
txt = open('proof_lead.md').read()
checks = {
 'P1 GG sign lemma5.1': 'generating function is `sum_{d odd} (...) z^d/d! = GG(z)/4`',
 'P2 Thm5.5 GG sign': '(u m_e part)/d! = (1/4)[z^d] GG(z)',
 'P3 Theorem 6.2': '**Theorem 6.2 (the law).**',
 'P3b Corollary 6.3': '**Corollary 6.3 (generating functions).**',
 'P4 branch v>0': '(the branch with `v > 0`;',
 'P5 (2) derivation': '(2) is the direct simplification',
 'P6 M_g strip': 'converges on `-1 < Re s < 2`',
 'P7 AC extension': 'so the identity extends to the\nwhole strip',
 'P8 sigma typo': '`sigma(p,q) = D_*` iff `(p,q) = D_*`)',
 'P9 4.2 opening': '**Proof.** Lemma 4.1 gives `K_d = -(1/2) int (L1-L2)^d m_e du',
 'P10 5.1 wording': 'The two sign flips cancel:',
 'P11 ranges 5.3': 'verified to 51-72 digits\nfor d = 1..11',
 'P12a notation Lemma 6.1': 'eps_{d-1} A_d/d!` (Lemma 6.1);',
 'P12b section 4 ref': 're-derived here (Section 4).',
 'P13 Cor(a) index': '= [u^w] (u^2 (1-sinh u)/cosh^2 u)` by Lemma 6.1.',
 'P14 d9 nuance': 'determined Newton fit on\n`k = 10..18` with zero spare points',
}
for k,v in checks.items():
    print(('OK ' if v in txt else 'MISSING '), k)
