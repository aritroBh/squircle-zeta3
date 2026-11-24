
import os, datetime
wd = "<project>/research/08_theory/SRF3_octahedron_radius_proof"
log = open(wd+"/progress.log","w")
entries = [
 "step_01: exact representation S=8(1-tau)^2 T(tau) (simplex + square + Laplace forms) verified vs H1 direct quadrature: 40 digits, 6 values of d [NUMERICAL]",
 "step_02: Goertzel/scan of normalized coeffs: single dominant frequency arg=1.787953 (edge), vertex arg 1.3536 absent [NUMERICAL]",
 "step_03: Pade-Borel poles noisy; superseded by annihilator/fits [NUMERICAL]",
 "step_04: edge-annihilator decays; beta=-1.5 sharply optimal => square-root singularity; vertex amplitude <0.2% [NUMERICAL]",
 "engine K=120: tau- and d-coefficients of S/S(1) to order 120 at ~460 dps (H1 Route-A engine, local copy) [NUMERICAL]",
 "fits K=120: C_tau=1.37661110 vs theory sqrt(3*pi*R/8)=1.37661108431311 (8 digits); phi=theta/2-pi to 1e-8; beta=-1.5000 [NUMERICAL]",
 "d-series: C_d=sqrt(2)*C_tau, phi_d=phi+3pi/2 (mixing factor e^{-z0bar}=sqrt(2)i) to 8 digits; mixing formula exact (1e-51) [NUMERICAL]",
 "theory: Borel singular part B(z) ~ D(1-z/z0)^{1/2}, D=(pi sqrt(3)/(2 sqrt(2))) sqrt(z0bar), from edge strips; s-integral contributes B(3/2,3/2)=pi/8 [PROVED in proof.md]",
]
for e in entries:
    log.write(f"[{datetime.date.today()}] {e}\n")
log.close()
print(open(wd+"/progress.log").read())
