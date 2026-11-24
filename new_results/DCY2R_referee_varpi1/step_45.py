import sys
sys.path.insert(0, '../DCY2_decay_first_correction')
from flint import arb, ctx
ctx.dps = 40
from dcy2_cell import *
from dcy2_sweep import window_cell, TR_cell, TL_cell, T0_cell

N = 1e9
Nlo, Nhi = N, N*1.0001
bracket_cell = window_cell(Nlo, Nhi)
tp_cell = TR_cell(Nlo,Nhi) + TL_cell(Nlo,Nhi)
y = y_of(N)
lhs_cellform = (2/PI)*y*S_sup(N)*(H_up_of(N)*bracket_cell + tp_cell)*arb(N)**2 + T0_cell(Nlo,Nhi)
print("pieces()/cell form at 1e9:", lhs_cellform.str(10))

lhs_pw = LHS_of(N)
print("LHS_of pointwise at 1e9:  ", lhs_pw.str(10))

# decompose LHS_of
br_pw = E_mom_of(N) + E_tailmom_of(N)
print("E_mom_of:", E_mom_of(N).str(6), " E_tailmom_of:", E_tailmom_of(N).str(6))
print("bracket_cell:", bracket_cell.str(6))
print("TR_cell:", TR_cell(Nlo,Nhi).str(6), " TL_cell:", TL_cell(Nlo,Nhi).str(6))
D = Delta_of(N)
print("eLamD*CH_R/rD:", (eLamD_up(N)*CH_R_of(N)/rDelta_of(N)).str(6))
print("LogM(-D)*CH_L/lamL:", (LogM_of(N,-D).exp()*(PI/2).sqrt()*2**(arb(3)/4)/lamL_of(N)).str(6))
