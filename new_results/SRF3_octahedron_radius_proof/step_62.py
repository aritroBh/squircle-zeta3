
import os, datetime
wd = "<project>/research/08_theory/SRF3_octahedron_radius_proof"
log = open(wd+"/progress.log","a")
log.write(f"[{datetime.date.today()}] verify_srf3.py final run: all PASS in 139s (V1 40.2 digits, V5 C rel err 1.4e-6@K60, V6 26.7 digits)\n")
log.write(f"[{datetime.date.today()}] SRF3 COMPLETE: proof.md (Thms A-E, lemma ledger), result.md, verify_srf3.py, K=120 artifacts\n")
log.close()
print(open(wd+"/progress.log").read())
print("\nFiles:")
for f in ["proof.md","result.md","verify_srf3.py","SRF3_hi_coeffs_K120.json","w_series_K120.json","engine_h1.py"]:
    print(" ", f, os.path.getsize(os.path.join(wd,f)), "bytes")
