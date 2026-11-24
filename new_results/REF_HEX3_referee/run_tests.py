import sys, time, sympy as sp, mpmath as mp
import direct, formula

DPS = int(sys.argv[1]) if (__name__ == "__main__" and len(sys.argv) > 1) else 60
mp.mp.dps = DPS

def frame(a, b, M):
    """corners (1,0),(a,b),(0,1) mapped by the 2x2 matrix M (det>0)."""
    pts = [(1, 0), (a, b), (0, 1)]
    return [(M[0][0] * x + M[0][1] * y, M[1][0] * x + M[1][1] * y) for x, y in pts]

S = sp.sqrt
HEX = {
  'H1_regular_w111': ([(S(3)/2/S(3), sp.Rational(1,2)/S(3)), (0, 1/S(3)), (-S(3)/2/S(3), sp.Rational(1,2)/S(3))], [1,1,1]),
  'H2_paper_example_w131': ([(1,0),(2,2),(0,1)], [1,3,1]),
  'H3_far_r2': (frame(2, sp.Rational(3,2), [[3,1],[0,sp.Rational(1,2)]]), [2,7,3]),
  'H4_r_half': (frame(sp.Rational(13,10), 2, [[1,-2],[1,1]]), [5,1,2]),
  'H5_irrational': (frame(S(2), S(3)-sp.Rational(1,2), [[1,S(5)],[0,1]]), [sp.Rational(3,2),4,1]),
  'H6_r4': (frame(4, sp.Rational(16,5), [[1,0],[0,1]]), [1,2,3]),
  'H7_r2_pow2logs': (frame(2, sp.Rational(3,2), [[2,1],[1,3]]), [1,5,3]),
}

def main():
  for name, (V, w) in HEX.items():
      eds = formula.edges_data(V, w)
      Vm = [(mp.mpf(sp.N(x, DPS + 20)), mp.mpf(sp.N(y, DPS + 20))) for x, y in V]
      wm = [mp.mpf(sp.N(x, DPS + 20)) for x in w]
      H = direct.setup(Vm, wm)
      t = time.time()
      tot_d = 0; tot_f = 0; tot_lit = 0; tot_0 = 0
      print(f'=== {name}: V={V} w={w}')
      for q in range(3):
          ed = eds[q]
          jet = direct.edge_jet_fast(H, q)
          Ef = formula.E_formula(ed, 'ell')
          El = formula.E_formula(ed, 'literal')
          tot_0 += jet[0]; tot_d += jet[3]; tot_f += Ef; tot_lit += El
          print(f'  edge {q}: r={ed["r"]} B={ed["B"]} L2={ed["L2"]}  direct E={mp.nstr(jet[3], 45)}')
          print(f'          formula-direct = {mp.nstr(Ef - jet[3], 5)}   literal-label formula - direct = {mp.nstr(El - jet[3], 5)}')
      Pex = 2*sum(mp.sqrt(ed['L2'] if isinstance(ed['L2'], mp.mpf) else mp.mpf(sp.N(ed['L2'], DPS+20))) for ed in eds)
      print(f'  P(1) direct - exact perimeter = {mp.nstr(2*tot_0 - Pex, 5)}')
      print(f"  P'''(1) direct = {mp.nstr(2*tot_d, 50)}")
      print(f"  P'''(1) table  = {mp.nstr(2*tot_f, 50)}   diff {mp.nstr(2*(tot_f-tot_d), 5)}   ({time.time()-t:.0f}s)")
      sys.stdout.flush()

if __name__=="__main__":
    main()
