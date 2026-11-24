import sys
sys.path.insert(0, '<project>/research/10_frontier/G1_zeta3_any_angle')
import mpmath as mp
import corner_series as cs

# convention test at c=1/2 (alpha = 2pi/3): stored COR5 value
stored = mp.mpf("0.010717362932289465418409576291224895855718975738464335760155612417841103619782418171396503742199126038436380620836337865304662784078999406937965739269009749615547401879802207222417678105864671523471037376024873555429679637766699958870900647844995948532283743592952407677723163070778394289744442351066844742744206712464681")
mp.mp.dps = 80
v_full = cs.corner_coeffs(2*mp.pi/3, 3)[3]
v_half = cs.corner_coeffs(mp.pi/3, 3)[3]
print("2pi/3 (full):", mp.nstr(v_full, 25), " err", mp.nstr(abs(v_full-stored),3))
print(" pi/3 (half):", mp.nstr(v_half, 25), " err", mp.nstr(abs(v_half-stored),3))
