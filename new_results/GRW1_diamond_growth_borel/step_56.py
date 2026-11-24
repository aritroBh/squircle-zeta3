import os, time
wd = "<project>/research/08_theory/GRW1_diamond_growth_borel"
log = open(wd+"/progress.log","a")
lines = [
 "step_02: identity I=int e^-u H K du verified (10 eta,q pairs, 30-44 digits); H_n exact rationals; H^[0] bounded on R+ (~1.4-1.5)",
 "step_03: H(z)=sqrt2 sqrt(1+e^-2z) confirmed (0/1e-81); H^[0](x)=sqrt2 sum c_k J_0(2sqrt(2kx)) confirmed to 1e-9",
 "step_04: M_{n,0}=n! eta(n) ✓; tail bound |M_{n,j}|=(n-j)!/(j!(j+1)^{n-j+1}) ✓; f_n=sum_j H_{n-j}M_{n,j} ✓ (3e-18..3e-39); tail/main=O(2^-n)",
 "step_05: growth law ratios f_n/law -> 1 all 4 classes (n<=180), Richardson const 1+-1.5e-7; f_n/(n!H_n)=1-2^-n...",
 "step_06/07: rhombus engine (operator formula, -log(1+eta) fix); f_1=-1 ✓; A,B swap symmetry ✓",
 "step_08-11: square 1/n term = (3-3pi(-1)^n)/(8n) identified (fits 1e-7); two-term Darboux law2 residual O(n^-2) const~4 (square) and rhombus (worst 4.8e-4 at safe n<=240)",
 "step_12/18: Pade-Borel [90/90] build (2.25e-358); direct PB digits 44.6/31.3/25.4/23.8 at 0.1/0.2/0.3/-0.2 (paper:44/31/25/30)",
 "step_20/21: R=B-H has radius 2rho=pi; H-split PB digits 57.0/43.7/35.4/42.7 at 0.1/0.2/0.3/-0.2; 34-57 digits over ±0.05..±0.3",
 "step_23: B(zeta)=sum Bj(zeta) integrals vs Pade vs series: agree 5e-18..2e-9 on real axis (zeta=±0.5..±3)",
 "step_25/26/27: Stokes jump J=2sqrt(pi)e^{-3i pi/4} eps^{1/2}e^{-pi/(2eps)}: ratios 0.9803,0.9889,0.9950,1.0003; phases ->0; R-jump ~ e^{-pi/eps} negligible; S_-(i eps)=I(i eps) to 1e-7..4e-12",
 "step_28: PB poles refined: all |z|>pi/2 (innermost 1.57175081); none near real axis",
 "step_29: S(eta)=I(eta) at ±0.5,±0.7,±0.9: 5.4-26.9 digits",
 "step_34: CERTIFIED (arb) I(eta) at 8 points, total radius <= 5.4e-62; match DE refs to 55+ digits",
]
for l in lines: log.write(l+"\n")
log.close()
print("logged", len(lines), "lines")
print(os.listdir(wd))
