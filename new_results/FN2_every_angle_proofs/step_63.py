import os
print("=== FINAL DELIVERABLES ===")
for f in ["proof.md","result.md","paper_proposal.tex","verify_fn2.py","verify_fn2.out",
          "progress.log","README_STEPS.md","residuals_main.json","residuals_K4.json","contact_results.json"]:
    print(f"  {f}: {os.path.getsize(f)} bytes")
print("\n=== proof.md section headers ===")
with open("proof.md") as f:
    for line in f:
        if line.startswith("## ") or line.startswith("**Status"):
            print(" ", line.rstrip())
