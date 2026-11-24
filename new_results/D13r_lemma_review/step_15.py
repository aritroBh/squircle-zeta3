work='research/08_theory/D13r_lemma_review'
with open(work+'/progress.log','a') as f:
    f.write("L1: HOLDS (hand proof: psi-series + |n+w|>=n+t termwise + integral test psi'<=1/t^2+1/t, |psi''|<=2/t^3+1/t^2; "
            "numeric 1659 pts @60dps: worst |g'|/G1b=0.80094308518 at x=9.85e7 (limit 0.8009431025), "
            "|g''|/G2b=0.62499997038 (limit 0.625), psi'/bound->1-, |psi''|/bound->1- as t->0, never exceeds). def of gp/gpp cross-checked vs finite diff of log g (resid 1.6e-27).\n")
print(open(work+'/progress.log').read())
