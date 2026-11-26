import Mathlib
/-!
# T12 -- edge law: parallelogram and affine-regular hexagon (prop:edgelaw)

Target id: T12.  Paper label: `prop:edgelaw` (main.tex lines 1671-1684; detailed proof in
`paper/ieee/proofs/sec_g4_edgelaw.tex`, Proposition 7 / Lemmas g4:lem:fl, g4:lem:laws).

Quoted claims (main.tex 1672-1683):

  "Let $e$ be an edge of a centrally symmetric convex polygon, from the corner $V_0$ to the
  corner $V_1$, of length $L$ and at distance $h$ from the center. The foot of the perpendicular
  from the center cuts the line of $e$ into signed lengths $L_0$ and $L_1=L-L_0$, measured from
  $V_0$ and $V_1$. If $\langle v_{j_0},x\rangle=0$ is the kink line through $V_0$ and
  $\langle v_{j_1},x\rangle=0$ the one through $V_1$, put $\xi_0=|\langle v_{j_0},V_1\rangle|$
  and $\xi_1=|\langle v_{j_1},V_0\rangle|$. Then the $\zeta(3)$ part of what the local terms of
  $e$ contribute to $P'''(1)$ is exactly $Z_e\,\zeta(3)$, where
  \[ Z_e=6\,\xi_0\xi_1\frac{h^2}{L}\Bigl(3\,\frac{\xi_0L_0+\xi_1L_1}{L}-2\Bigr). \]
  ... With $\mathcal E=\sum_eh_e^2/L_e$, every parallelogram has
  \[ P(1+d)=P(1)+\frac{\pi^2\mathcal E}{6}d^2+\mathcal E\Bigl(\zeta(3)-\frac{\pi^2}{3}\Bigr)d^3+O(d^4), \]
  and for every affine-regular hexagon $P(1+d)$ is
  \[ 2^{\frac{d}{1+d}}\Bigl[P(1)+\frac{\pi^2\mathcal E}{24}d^2-\frac{\mathcal E}{8}
     \Bigl(\zeta(3)+\frac{2\pi^2}{3}\Bigr)d^3\Bigr]+O(d^4). \]"

What is formalized: the algebraic consistency step between the edge law and the two displays.
The polygon is $K_1=\{\sum_j|\langle v_j,x\rangle|\le1\}$ in $\mathbb R^2$ (modelled as `ℝ × ℝ`
with the Euclidean dot product). All edge data are computed from the corners:
* $L=|V_1-V_0|$;
* $L_0=-\langle V_0,V_1-V_0\rangle/L$ (signed length from $V_0$ to the foot $F$ of the
  perpendicular from the origin, positive towards $V_1$; this matches main.tex line 1695,
  $q_0=V_0\cdot(V_1-V_0)=-LL_0$), $L_1=L-L_0$;
* $F=V_0+(L_0/L)(V_1-V_0)$ and $h=|F|$;
* $\xi_0=|\langle v_{j_0},V_1\rangle|$, $\xi_1=|\langle v_{j_1},V_0\rangle|$ where $v_{j_0}$
  ($v_{j_1}$) is the vector whose kink line passes through $V_0$ ($V_1$).
The corners are not chosen: they are any points satisfying the paper's defining conditions
(on the boundary $\sum_j|\langle v_j,x\rangle|=1$ and on the named kink line).

* Parallelogram: $m=2$, $v_1,v_2$ non-parallel; $C$ is a boundary point on the kink line
  $\langle v_1,x\rangle=0$ and $D$ one on $\langle v_2,x\rangle=0$; the four edges are
  $C\to D\to -C\to -D\to C$.
* Affine-regular hexagon: $m=3$; corners $C_0,C_1,C_2,-C_0,-C_1,-C_2$ with $C_0,C_2$ non-parallel
  and $C_1=C_0+C_2$ (equivalently, a linear image of the regular hexagon, whose corners satisfy
  $R_1=R_0+R_2$; this is the characterization of Lemma g4:lem:fl(c) of the supplement); $C_i$ lies
  on the kink line $\langle v_i,x\rangle=0$ and on the boundary of $K_1$.

The $\zeta(3)$-coefficient of $P'''(1)$ implied by each display is extracted by treating
$\zeta(3)$ as a formal real parameter `z3` (the displays are affine in it): the $\zeta(3)$-part of
$P'''(1)=3!\,[d^3]$ is `iteratedDeriv 3 (disp z3 - disp 0) 0`, which must equal
$(\sum_eZ_e)\,z3$.  Here `z3` is a universally quantified real (formal coefficient identity),
not `riemannZeta 3`.  The factor $2^{d/(1+d)}$ is kept in the hexagon display (`Real.rpow`; the
value of `d/(1+d)` at `d = -1` is irrelevant since only the germ at `d = 0` is used).

No hypothesis is strengthened beyond the paper: "pairwise non-parallel $v_j$" is the paper's
setting; non-parallel $C_0,C_2$ is part of being a (non-degenerate) hexagon.

NOT COVERED:
* The analytic content of the two displays (that the perimeter $P(p)$ of
  $K_p=\{\sum_j|\langle v_j,x\rangle|^p\le1\}$ really has these expansions with an $O(d^4)$
  error).  Mathlib has no usable perimeter/arc-length for $\partial K_p$; one could use
  `μH[1] (frontier K_p)`, but proving it is far beyond reach, so it is left out.
* The general claim that the $\zeta(3)$ part of the local terms of an arbitrary edge equals
  $Z_e\zeta(3)$ (needs the notions "local terms" and "formal $\zeta(3)$ part" of a Beta-moment
  expansion of the length element), and "if the polygon is fully local this is all of its
  $\zeta(3)$".
* The remarks after the proposition ($\sum_eZ_e=C_m$, the regular-hexagon value of $P'''(1)$).
-/
-- BEGIN STATEMENTS
namespace T12

/-- Euclidean inner product on `ℝ × ℝ`. -/
def dot (x y : ℝ × ℝ) : ℝ := x.1 * y.1 + x.2 * y.2

/-- Determinant of two plane vectors (nonzero iff they are not parallel). -/
def det (x y : ℝ × ℝ) : ℝ := x.1 * y.2 - x.2 * y.1

/-- Edge length $L=|V_1-V_0|$. -/
noncomputable def edgeLen (V0 V1 : ℝ × ℝ) : ℝ := Real.sqrt (dot (V1 - V0) (V1 - V0))

/-- Signed length $L_0$ from $V_0$ to the foot of the perpendicular from the origin. -/
noncomputable def edgeL0 (V0 V1 : ℝ × ℝ) : ℝ := -dot V0 (V1 - V0) / edgeLen V0 V1

/-- $L_1=L-L_0$. -/
noncomputable def edgeL1 (V0 V1 : ℝ × ℝ) : ℝ := edgeLen V0 V1 - edgeL0 V0 V1

/-- Foot of the perpendicular from the origin to the line of the edge. -/
noncomputable def foot (V0 V1 : ℝ × ℝ) : ℝ × ℝ :=
  V0 + (edgeL0 V0 V1 / edgeLen V0 V1) • (V1 - V0)

/-- Distance $h$ of the edge line from the center. -/
noncomputable def edgeH (V0 V1 : ℝ × ℝ) : ℝ := Real.sqrt (dot (foot V0 V1) (foot V0 V1))

/-- An oriented edge `V0 → V1` together with the vectors `u0 = v_{j0}`, `u1 = v_{j1}` whose kink
lines pass through `V0` and `V1`. -/
structure Edge where
  u0 : ℝ × ℝ
  u1 : ℝ × ℝ
  V0 : ℝ × ℝ
  V1 : ℝ × ℝ

/-- The paper's $Z_e=6\xi_0\xi_1\frac{h^2}{L}\bigl(3\frac{\xi_0L_0+\xi_1L_1}{L}-2\bigr)$ with
$\xi_0=|\langle v_{j_0},V_1\rangle|$, $\xi_1=|\langle v_{j_1},V_0\rangle|$. -/
noncomputable def Ze (e : Edge) : ℝ :=
  let ξ0 := |dot e.u0 e.V1|
  let ξ1 := |dot e.u1 e.V0|
  let L := edgeLen e.V0 e.V1
  let h := edgeH e.V0 e.V1
  6 * ξ0 * ξ1 * (h ^ 2 / L) *
    (3 * (ξ0 * edgeL0 e.V0 e.V1 + ξ1 * edgeL1 e.V0 e.V1) / L - 2)

/-- The summand $h_e^2/L_e$ of $\mathcal E$. -/
noncomputable def hsqOverL (e : Edge) : ℝ := edgeH e.V0 e.V1 ^ 2 / edgeLen e.V0 e.V1

/-- $\sum_e Z_e$ over a list of edges. -/
noncomputable def sumZ (es : List Edge) : ℝ := (es.map Ze).sum

/-- $\mathcal E=\sum_e h_e^2/L_e$ over a list of edges. -/
noncomputable def sumE (es : List Edge) : ℝ := (es.map hsqOverL).sum

/-- Corner data of the parallelogram $K_1=\{|\langle v_1,x\rangle|+|\langle v_2,x\rangle|\le1\}$:
`C` is the boundary point on the kink line $\langle v_1,x\rangle=0$, `D` the one on
$\langle v_2,x\rangle=0$. -/
def IsParallelogram (v1 v2 C D : ℝ × ℝ) : Prop :=
  det v1 v2 ≠ 0 ∧
  dot v1 C = 0 ∧ |dot v1 C| + |dot v2 C| = 1 ∧
  dot v2 D = 0 ∧ |dot v1 D| + |dot v2 D| = 1

/-- The four edges $C\to D\to-C\to-D\to C$, each with the kink vectors through its endpoints. -/
def parEdges (v1 v2 C D : ℝ × ℝ) : List Edge :=
  [⟨v1, v2, C, D⟩, ⟨v2, v1, D, -C⟩, ⟨v1, v2, -C, -D⟩, ⟨v2, v1, -D, C⟩]

/-- Corner data of an affine-regular hexagon
$K_1=\{\sum_{i=0}^2|\langle v_i,x\rangle|\le1\}$ with corners $C_0$, $C_1=C_0+C_2$, $C_2$ and
their negatives; $C_i$ lies on the kink line $\langle v_i,x\rangle=0$ and on $\partial K_1$. -/
def IsAffRegHexagon (v0 v1 v2 C0 C2 : ℝ × ℝ) : Prop :=
  det C0 C2 ≠ 0 ∧ det v0 v1 ≠ 0 ∧ det v1 v2 ≠ 0 ∧ det v0 v2 ≠ 0 ∧
  dot v0 C0 = 0 ∧ dot v1 (C0 + C2) = 0 ∧ dot v2 C2 = 0 ∧
  |dot v0 C0| + |dot v1 C0| + |dot v2 C0| = 1 ∧
  |dot v0 (C0 + C2)| + |dot v1 (C0 + C2)| + |dot v2 (C0 + C2)| = 1 ∧
  |dot v0 C2| + |dot v1 C2| + |dot v2 C2| = 1

/-- The six edges $C_0\to C_1\to C_2\to-C_0\to-C_1\to-C_2\to C_0$ ($C_1=C_0+C_2$). -/
def hexEdges (v0 v1 v2 C0 C2 : ℝ × ℝ) : List Edge :=
  let C1 := C0 + C2
  [⟨v0, v1, C0, C1⟩, ⟨v1, v2, C1, C2⟩, ⟨v2, v0, C2, -C0⟩,
   ⟨v0, v1, -C0, -C1⟩, ⟨v1, v2, -C1, -C2⟩, ⟨v2, v0, -C2, C0⟩]

/-- Parallelogram display with $\zeta(3)$ replaced by the formal parameter `z3`:
$P+\frac{\pi^2\mathcal E}6d^2+\mathcal E(z_3-\frac{\pi^2}3)d^3$. -/
noncomputable def parDisp (P E z3 d : ℝ) : ℝ :=
  P + Real.pi ^ 2 * E / 6 * d ^ 2 + E * (z3 - Real.pi ^ 2 / 3) * d ^ 3

/-- Hexagon display with $\zeta(3)$ replaced by `z3`:
$2^{d/(1+d)}\bigl[P+\frac{\pi^2\mathcal E}{24}d^2-\frac{\mathcal E}8(z_3+\frac{2\pi^2}3)d^3\bigr]$. -/
noncomputable def hexDisp (P E z3 d : ℝ) : ℝ :=
  (2 : ℝ) ^ (d / (1 + d)) *
    (P + Real.pi ^ 2 * E / 24 * d ^ 2 - E / 8 * (z3 + 2 * Real.pi ^ 2 / 3) * d ^ 3)

/-- Parallelogram, edge law summed: with edge data computed from the corners,
$\sum_eZ_e=6\,\mathcal E$ (= $3!$ times the $\zeta(3)$-coefficient $\mathcal E$ of $d^3$). -/
def stmt_1 : Prop :=
  ∀ v1 v2 C D : ℝ × ℝ, IsParallelogram v1 v2 C D →
    sumZ (parEdges v1 v2 C D) = 6 * sumE (parEdges v1 v2 C D)

/-- Parallelogram: the $\zeta(3)$-part of $P'''(1)=3!\,[d^3]$ read off the paper's display equals
$(\sum_eZ_e)\,\zeta(3)$ (formal `z3`). -/
def stmt_2 : Prop :=
  ∀ v1 v2 C D : ℝ × ℝ, IsParallelogram v1 v2 C D → ∀ P z3 : ℝ,
    iteratedDeriv 3
        (fun d => parDisp P (sumE (parEdges v1 v2 C D)) z3 d
                  - parDisp P (sumE (parEdges v1 v2 C D)) 0 d) 0
      = sumZ (parEdges v1 v2 C D) * z3

/-- Affine-regular hexagon, edge law summed: $\sum_eZ_e=6\cdot(-\mathcal E/8)$. -/
def stmt_3 : Prop :=
  ∀ v0 v1 v2 C0 C2 : ℝ × ℝ, IsAffRegHexagon v0 v1 v2 C0 C2 →
    sumZ (hexEdges v0 v1 v2 C0 C2) = 6 * (-(sumE (hexEdges v0 v1 v2 C0 C2) / 8))

/-- Affine-regular hexagon: the $\zeta(3)$-part of $3!\,[d^3]$ of the paper's display (with the
factor $2^{d/(1+d)}$) equals $(\sum_eZ_e)\,\zeta(3)$ (formal `z3`). -/
def stmt_4 : Prop :=
  ∀ v0 v1 v2 C0 C2 : ℝ × ℝ, IsAffRegHexagon v0 v1 v2 C0 C2 → ∀ P z3 : ℝ,
    iteratedDeriv 3
        (fun d => hexDisp P (sumE (hexEdges v0 v1 v2 C0 C2)) z3 d
                  - hexDisp P (sumE (hexEdges v0 v1 v2 C0 C2)) 0 d) 0
      = sumZ (hexEdges v0 v1 v2 C0 C2) * z3

end T12
-- END STATEMENTS

open T12 in
lemma T12.Ze_eq_of_xi (e : Edge) (c : ℝ) (h0 : |dot e.u0 e.V1| = c) (h1 : |dot e.u1 e.V0| = c) :
    Ze e = 6 * (c ^ 2 * (3 * c - 2)) * hsqOverL e := by
  unfold Ze hsqOverL
  simp only [h0, h1]
  by_cases hL : edgeLen e.V0 e.V1 = 0
  · simp [hL]
  · have : edgeL0 e.V0 e.V1 + edgeL1 e.V0 e.V1 = edgeLen e.V0 e.V1 := by
      unfold edgeL1; ring
    have h2 : c * edgeL0 e.V0 e.V1 + c * edgeL1 e.V0 e.V1 = c * edgeLen e.V0 e.V1 := by
      rw [← mul_add, this]
    rw [h2]
    field_simp

open T12 in
lemma T12.sumZ_eq (es : List Edge) (k : ℝ) (h : ∀ e ∈ es, Ze e = k * hsqOverL e) :
    sumZ es = k * sumE es := by
  unfold sumZ sumE
  induction es with
  | nil => simp
  | cons a t ih =>
    simp only [List.map_cons, List.sum_cons]
    rw [h a (by simp), ih (fun e he => h e (by simp [he]))]
    ring

lemma T12.iterDeriv_aux (g : ℝ → ℝ) (hg : ContDiffAt ℝ 3 g 0) :
    iteratedDeriv 3 (fun d => g d * d ^ 3) 0 = 6 * g 0 := by
  have h2 : ContDiffAt ℝ 3 (fun d : ℝ => d ^ 3) 0 := by fun_prop
  rw [iteratedDeriv_fun_mul hg h2]
  simp [Finset.sum_range_succ]
  norm_num [Nat.factorial]
  ring

open T12 in
lemma T12.par_edges (v1 v2 C D : ℝ × ℝ) (h : IsParallelogram v1 v2 C D) :
    sumZ (parEdges v1 v2 C D) = 6 * sumE (parEdges v1 v2 C D) := by
  obtain ⟨-, h1, h2, h3, h4⟩ := h
  have e1 : |dot v2 C| = 1 := by rw [h1] at h2; simpa using h2
  have e2 : |dot v1 D| = 1 := by rw [h3] at h4; simpa using h4
  have hd : ∀ x y : ℝ × ℝ, dot x (-y) = - dot x y := by
    intro x y; simp [dot]; ring
  have hd2 : ∀ x y : ℝ × ℝ, dot (-x) y = - dot x y := by
    intro x y; simp [dot]; ring
  have := T12.sumZ_eq (parEdges v1 v2 C D) 6 ?_
  · exact this
  intro e he
  simp only [parEdges, List.mem_cons, List.not_mem_nil, or_false] at he
  rcases he with rfl | rfl | rfl | rfl
  all_goals
    rw [T12.Ze_eq_of_xi _ 1 (by simp [hd, hd2, e1, e2]) (by simp [hd, hd2, e1, e2])]
    norm_num

theorem T12.proof_1 : T12.stmt_1 := by
  intro v1 v2 C D h
  exact T12.par_edges v1 v2 C D h

theorem T12.proof_2 : T12.stmt_2 := by
  intro v1 v2 C D h P z3
  have h1 := T12.par_edges v1 v2 C D h
  have : (fun d => T12.parDisp P (T12.sumE (T12.parEdges v1 v2 C D)) z3 d
                  - T12.parDisp P (T12.sumE (T12.parEdges v1 v2 C D)) 0 d)
      = fun d => (T12.sumE (T12.parEdges v1 v2 C D) * z3) * d ^ 3 := by
    funext d; unfold T12.parDisp; ring
  rw [this, iteratedDeriv_const_mul_field]
  simp [iteratedDeriv_fun_pow_zero, h1, Nat.factorial]
  ring

open T12 in
lemma T12.hex_edges (v0 v1 v2 C0 C2 : ℝ × ℝ) (h : IsAffRegHexagon v0 v1 v2 C0 C2) :
    sumZ (hexEdges v0 v1 v2 C0 C2) = 6 * (-(sumE (hexEdges v0 v1 v2 C0 C2) / 8)) := by
  obtain ⟨-, -, -, -, a0, a1, a2, b0, b1, b2⟩ := h
  have add : ∀ x y z : ℝ × ℝ, dot x (y + z) = dot x y + dot x z := by
    intro x y z; simp [dot]; ring
  have hd : ∀ x y : ℝ × ℝ, dot x (-y) = - dot x y := by
    intro x y; simp [dot]; ring
  have hd2 : ∀ x y : ℝ × ℝ, dot (-x) y = - dot x y := by
    intro x y; simp [dot]; ring
  simp only [add] at a1 b1
  rw [a0] at b0
  rw [a0, a1, a2] at b1
  rw [a2] at b2
  simp only [abs_zero, zero_add, add_zero] at b0 b1 b2
  have e1 : dot v1 C2 = - dot v1 C0 := by linarith
  rw [e1, abs_neg] at b2
  have hp : |dot v0 C2| = 1 / 2 := by linarith
  have hq : |dot v1 C0| = 1 / 2 := by linarith
  have hr : |dot v2 C0| = 1 / 2 := by linarith
  have := T12.sumZ_eq (hexEdges v0 v1 v2 C0 C2) (6 * (-(1 / 8))) ?_
  · rw [this]; ring
  intro e he
  simp only [hexEdges, List.mem_cons, List.not_mem_nil, or_false] at he
  rcases he with rfl | rfl | rfl | rfl | rfl | rfl
  all_goals
    rw [T12.Ze_eq_of_xi _ (1 / 2) (by simp [add, hd, hd2, a0, a2, e1, abs_neg, hp, hq, hr])
      (by simp [add, hd, hd2, a0, a2, e1, abs_neg, hp, hq, hr])]
    norm_num

theorem T12.proof_3 : T12.stmt_3 := by
  intro v0 v1 v2 C0 C2 h
  exact T12.hex_edges v0 v1 v2 C0 C2 h

theorem T12.proof_4 : T12.stmt_4 := by
  intro v0 v1 v2 C0 C2 h P z3
  have h1 := T12.hex_edges v0 v1 v2 C0 C2 h
  set E := T12.sumE (T12.hexEdges v0 v1 v2 C0 C2)
  have : (fun d => T12.hexDisp P E z3 d - T12.hexDisp P E 0 d)
      = fun d => (2 : ℝ) ^ (d / (1 + d)) * d ^ 3 * (-(E / 8) * z3) := by
    funext d; unfold T12.hexDisp; ring
  have hg : ContDiffAt ℝ 3 (fun d : ℝ => (2 : ℝ) ^ (d / (1 + d))) 0 := by
    have : (fun d : ℝ => (2 : ℝ) ^ (d / (1 + d))) = fun d => Real.exp (Real.log 2 * (d / (1 + d))) := by
      funext d; rw [Real.rpow_def_of_pos (by norm_num)]
    rw [this]
    refine ContDiffAt.exp (ContDiffAt.mul contDiffAt_const (ContDiffAt.div contDiffAt_id (contDiffAt_const.add contDiffAt_id) ?_))
    simp
  rw [this, iteratedDeriv_mul_const_field, h1]
  have := T12.iterDeriv_aux _ hg
  rw [this]
  simp
  ring
