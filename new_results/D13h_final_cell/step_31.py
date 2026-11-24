r = open('result.md').read()
r = r.replace(
 "| final cell C (J=192) | 49.56259362936965654550929024126328268240614812466812154794657124742081 | 2.7e-69 | 68 |\n"
 "| final cell C (J=48) | 53.61568583421857309778589251677539452160470915959575395741864585441471 | 2.8e-69 | 68 |",
 "| final cell C (J=192) | 49.56259792396024857830698314833209853198249410057901442521712151382263 | 2.0e-69 | 68 |\n"
 "| final cell C (J=96) | 50.85487159316713887613912557631075219361998757785207777044262300667321 | 4.5e-69 | 68 |\n"
 "| final cell C (J=48) | 53.61569012880916513058358542384421037118105513550664683468919612081652 | 2.6e-69 | 68 |")
r = r.replace("`step_10.py` (J=48, C=53.6157) and `step_10_J192.py` (J=192, **C=49.5626**)",
              "`step_10.py` (J=48, C=53.61569), `step_10_J96.py` (J=96, C=50.85487), `step_10_J192.py` (J=192, **C=49.56260**)")
r = r.replace("| D | 11.4423 | 11.4423 | step_38: ~2.44 (ours uses crude |E| bound) |",
              "| D | 11.44228 | 11.44228 | step_38: ~2.44 (ours uses crude |E| bound) |")
r = r.replace("| **C** | **53.6157** | **49.5626** | pointwise limit C(1e9) = 51.16 decreasing |",
              "| **C** | **53.61569** | **49.56260** | pointwise limit C(1e9) = 51.16 decreasing |")
r = r.replace("  $PY step_10_J192.py   # FINAL CELL J=192: C = 49.5626 for all N >= 1e10         (~2 s)",
              "  $PY step_10_J96.py    # FINAL CELL J=96 : C = 50.8549 for all N >= 1e10         (~1 s)\n  $PY step_10_J192.py   # FINAL CELL J=192: C = 49.5626 for all N >= 1e10         (~2 s)")
r = r.replace("  $PY step_10.py        # FINAL CELL J=48 : C = 53.6157 for all N >= 1e10         (<1 s)",
              "  $PY step_10.py        # FINAL CELL J=48 : C = 53.61569 for all N >= 1e10        (<1 s)")
open('result.md','w').write(r)
print("result.md updated")
with open('progress.log','a') as f:
    f.write("turn11: result.md updated with final post-d5-fix enclosures. DONE.\n")
