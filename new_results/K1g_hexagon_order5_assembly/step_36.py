import shutil, os
k1f = '<project>/research/08_theory/K1f_hexagon_order5'
WD = '<project>/research/08_theory/K1g_hexagon_order5_assembly'
# copy step_31.py (atom evaluator) and esum5.py into workspace for reuse
for f in ['step_31.py','esum5.py','atom_table.pkl']:
    shutil.copy(os.path.join(k1f,f), os.path.join(WD,'atom_eval.py' if f=='step_31.py' else f))
print("copied atom evaluator + esum5 + atom_table")
# quick test of atom evaluator
import importlib.util, sys
spec = importlib.util.spec_from_file_location("atom_eval", os.path.join(WD,'atom_eval.py'))
ae = importlib.util.module_from_spec(spec); spec.loader.exec_module(ae)
import mpmath as mp
mp.mp.dps=80
print("H1(-1) =", mp.nstr(ae.atom_eval('H1',-1),25), " expected -log(2) =", mp.nstr(-mp.log(2),25))
print("Li5(-1) =", mp.nstr(ae.atom_eval('Li5',-1),25))
print("  check -15/16 z5 =", mp.nstr(-mp.mpf(15)/16*mp.zeta(5),25))
print("Li2(-1) =", mp.nstr(ae.atom_eval('Li2',-1),25), " expected -z2/2 =", mp.nstr(-mp.zeta(2)/2,25))
