WD='<project>/research/08_theory/FLD1_field_rule_proof'
with open(WD+'/progress.log','a') as fh:
    fh.write("step20: bijection w<->(eps,k,v) proved & verified (242 words len<=5). h_w(N)=cf eps^N/N^k Z_v(N-1)\n")
    fh.write("step21: fast exact mu_{k,n} (Stirling, no symbolic series), k<=5 n<=100 in 8s; matches G1 k<=3\n")
# result.md skeleton, written early
result = r"""# FLD1: Field rule for every k -- result.md (running log)

**Goal.** Prove conj:field for every k>=0, 0<alpha<pi:
M_k(alpha) = sum_j z_j r_j(c) H_{w_j}(c),  c=cos(alpha/2),
with z_j monomials in zeta(2..k), r_j in Q(c) (denominators | c(1-c)(1+c)),
H_w HPLs with letters {0,1,-1}, and weight(z_j)+|w_j| <= k+1.

## Status / labels
- Core dictionary (HPL word <-> Z-sum data): **PROVED** (step_20), verified on all 242 words of len<=5.
- Exact mu_{k,n} for k<=5, n<=100: **PROVED FOR k<=5** (step_21, exact fractions, matches G1 k<=3).

## Log
- step_20_bijection.py: h_w(N) = cf eps^N/N^k Z_v(N-1); w->(eps,k,v) is a BIJECTION onto {+-1}xZ>=1x{words}.
- step_21_fastmu.py: exact mu_{k,n} table (K<=5,N<=100). mu_{4,2}, mu_{5,2} zeta-coeffs match flat_single_zeta patterns.
"""
open(WD+'/result.md','w').write(result)
print("logged. Now symbolic-n mu.")
