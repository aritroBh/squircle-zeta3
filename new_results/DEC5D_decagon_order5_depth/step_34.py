
import os
cwd = os.getcwd()
rmd = """# DEC5D — the regular decagon's P⁽⁵⁾(1): depth two needed? — status report

## Bottom line (to be filled)

Work in progress. This file is updated as steps complete.

## Plan / evidence table

| # | Check | Result | Residual / agreement | Digits | Producer |
|---|---|---|---:|---:|---|
| 1 | P⁽⁵⁾(1) two quadrature levels | pending | | | dec_engine.py |
"""
open(cwd+"/result.md","w").write(rmd)
print("result.md skeleton written")
# quick correctness run at low precision
import importlib, sys
sys.path.insert(0, cwd)
import dec_engine
res = dec_engine.run(40, "p5_run40.pkl")
