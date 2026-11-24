import time
def log_progress(line):
    with open('progress.log','a') as f: f.write(time.strftime('%H:%M:%S ')+line+'\n')
log_progress("BOUNDARY STRUCTURE (EXACT, Fraction engine): On stratum {a1=0,a2=a4,a5=a3}: C4=Q7'=Q11=0 IDENTICALLY, and ALL W12 irreducibles (z2^2z35,z2z37,CLASS,z2118) vanish EXACTLY (tested t=1,2,3,a3=1,2). Same for reversal stratum {a5=0,a2=a4,a1=a3}. Also faces {a2=0},{a3=0},{a4=0}: ALL irreducibles (W8,W10,W11,W12) vanish at all tested pts (a2a3a4 factor divides every irreducible coordinate). {a1=0}/{a5=0} OFF-stratum: irreducibles present. => Irreducible-free-through-W12 locus in CLOSED orthant = faces a2,a3,a4=0 PLUS strata {a1=0,a2=a4,a5=a3},{a5=0,a2=a4,a1=a3}; in OPEN orthant it is EMPTY.")
print(open('progress.log').read())
