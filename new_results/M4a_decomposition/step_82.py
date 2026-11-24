import os
want = ["ibp_reduce_k4.out","omega4.json","quad_90.json","m4_neville.json","gpl_w5.tsv","gpl_w4pp.tsv",
        "gpl_w5_check.json","gpl_w4pp_check.json","wv_w5.pkl","wv_w4pp.pkl","wv_w5_check.json",
        "wv_w4pp_check.json","hp_results.json","W4pp_800.txt","W5_600.txt","M4_800.txt","W4pp_1060.txt",
        "pslq_w4pp_700.json","pslq2_w4pp_700.json","pslq2_w4pp_1000.json","w4pp_closed_form.json",
        "result.md","progress.log","step_quad_pieces.py","step_neville_m4.py","gexpr.py","g_eval.py",
        "step_wordvectors.py","step_wv_check.py","step_wv_check5.py","step_hp_probe.py","step_hp.py",
        "step_pslq.py","step_pslq2.py","step_assemble_w4pp.py"]
for f in want:
    print(("OK " if os.path.exists(f) else "MISSING ") + f, os.path.getsize(f) if os.path.exists(f) else "")
