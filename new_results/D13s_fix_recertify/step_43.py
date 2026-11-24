
import mpmath as mp
mp.mp.dps = 30
pi = mp.pi
# From verify_TM_final.py output: TM_exact/sc etc (raw), certified C-parts (with 2/pi already)
data = {
 1e8:  dict(TM=0.250308597541, TR=0.228539735347, TLm=0.130082549743),
 1e10: dict(TM=0.226388652821, TR=0.177994533723, TLm=0.142694230242),
 1e12: dict(TM=0.208096065221, TR=0.152420799638, TLm=0.131444629283),
 1e14: dict(TM=0.193606409275, TR=0.137862737224, TLm=0.136051013855),
}
bnd = dict(TM=0.153370965962, TR=0.351962709577, TLm=1.92383977815)
print("piece: certified C-part vs (2/pi)*truth/sc  (margins)")
for N, d in data.items():
    line = f"N={N:g}: "
    for k in ['TM','TR','TLm']:
        truth_cpart = (2/pi)*d[k]
        line += f"{k} {mp.nstr(bnd[k],8)} > {mp.nstr(truth_cpart,8)} (x{mp.nstr(bnd[k]/truth_cpart,6)})   "
    print(line)
# specifically TM at N=1e10: the tightest
print()
print("TM at N=1e10: certified", mp.nstr(bnd['TM'],12), "vs truth (2/pi)*0.226388652821 =", mp.nstr((2/pi)*mp.mpf('0.226388652821'),12),
      " margin", mp.nstr(bnd['TM']/((2/pi)*mp.mpf('0.226388652821')),8))
with open('research/08_theory/D13s_fix_recertify/D13s_work/progress.log','a') as f:
    f.write("D13s turn10: FIXED final cell: C=47.07012493 (buggy 49.56259792; review corrected-with-old-breg 49.94909).\n"
            " TMreg 0.12237233 (was 0.10290245, x2^{1/4}); TM C-part 0.15337 > truth 0.14412 @1e10 (margin 1.064) CONFIRMED;\n"
            " TR 0.35196>0.11339 (3.1x); TLm 1.92384>0.09082 (21x). Btrue box [0.16666+/-8e-6], B_abs=0.16666814 <= breg_abs=0.75 (old harmless);\n"
            " d6 0.37024 (was 1.666), Dreg 8.56332 (was 11.44228). TLfN (pi/2)^1.5 restored: 5.5e-19229. TL0N proved form: 1.4e-2666376434.\n")
