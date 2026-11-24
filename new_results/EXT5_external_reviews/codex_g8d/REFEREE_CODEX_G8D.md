# Independent targeted second check: G8B′ and W₃

**Overall verdict: CORRECT AFTER REPAIR.** I find no counterexample to Corollary G8B′ and no failure of its pullback obstruction. Its hypotheses suffice, including for an arc of the unit circle. The revised W₃ regularisation and all requested counts are correct. Before publication, repair the inaccurate algebraic/genus scope sentences and two small proof-wording errors identified below. These repairs require no new hypothesis or mathematical argument for the stated corollary.

All line numbers refer to the supplied files. I did not edit them. My independent implementation is `check_g8d.py`; its successful output is `check_g8d.out`. It imports no supplied verification code. Reproduce with:

```sh
'<project>/.venv/bin/python' check_g8d.py
```

## 1. Exact statement and the angle case — CORRECT

**Locations:** `sec_g8_M3func.tex:239–244`, `:248–250`, `:277`.

The assumptions in line 242 are sufficient. At the accumulation point, φ is holomorphic and φ(t*) lies strictly inside D\{0}. Thus a disc B about t* has φ(B)⊂D\{0}. The chosen finite presentation of C excludes only finitely many coefficient poles; deleting those points from Γ leaves t* an accumulation point. A meromorphic identity on B follows. One can subsequently choose a smaller disc avoiding the entire finite alphabet. There is no need for Γ to contain an interval, for t* to belong to Γ, or for φ′(t*) to be nonzero.

For t=e^{iα/2}, any interior α* of J gives t* in the arc and in U. Here t*≠0, the only finite pole of φ(t)=(t+t^{-1})/2 is 0, and φ(t*)=cos(α*/2)∈(0,1). Every point of the arc is an accumulation point. All hypotheses hold. The claimed variants t=e^{iα/(2N)} and t=tan(α/4) also satisfy them with the displayed rational φ.

This conclusion is about equality on an angle interval with the specified branch presentations, not about a value at one angle.

## 2. Proof, step by step — CORRECT AFTER REPAIR

### (i) Meromorphic identity theorem — CORRECT AFTER REPAIR

**Location:** `sec_g8_M3func.tex:250`.

The meromorphic identity argument works even when t* is initially a pole of a presented coefficient. If F=C−𝔐∘φ has a pole of order m at t*, then (t−t*)^mF is holomorphic and has zeros accumulating at t*; the identity theorem makes it zero. Thus a nonzero meromorphic function cannot have such an accumulation of zeros, including at a pole.

There is a small mismatch between the hypothesis and the proof wording: the hypothesis excludes poles of **any coefficient**, whereas Step 1 says “not a pole of C”. Cancellation of coefficient poles can make those sets different. Use the set actually supplied by the hypothesis.

**Exact repair:** replace the sentence beginning “For t∈Γ∩B” through the next assertion of accumulation by:

> For t∈Γ∩B that is not a pole of any coefficient h_i of the chosen presentation of C, we have φ(t)∈(0,1), so 𝔐(φ(t))=M₃(φ(t))=C(t). The excluded coefficient poles form a finite set, so these equality points still accumulate at t*. Thus the meromorphic function C−𝔐∘φ on B has zeros accumulating at t*.

The remainder of Step 1 then applies exactly as written. This is a repair of the justification, not of the statement.

### (ii) Pullback of H_w and its top component — CORRECT

**Locations:** `sec_g8_M3func.tex:114–120`, `:144–155`, `:252–265`.

The chain rule in line 258 has the correct signs:

- H₀ differentiates against dlog φ;
- H₁ differentiates against **−**dlog(φ−1);
- H₋₁ differentiates against dlog(φ+1).

Consequently the sign is s_w=(−1)^{number of 1's in w}, and reversal of the word puts the outermost letter last. The induction uses only regular suffixes, which are regular. Its empty-word case is valid. Integration constants are length zero; changing the base point changes only shorter components. No global inverse to φ, unramifiedness assumption, or transport of a tangential base point through φ is being used.

There is no missing independent infinity letter. A rational function modulo nonzero constants is determined by its finite orders, and its order at infinity is minus their sum. In the angle substitution, all three pulled-back letters have a finite pole of order one at t=0 and a pole of order one at infinity. The former contributes **−e₀** to every letter; the latter is encoded by the sum of finite orders. Locally one avoids t=0; globally its letter must remain in the alphabet, as it does here. Constant factors such as −1 and 1/2 affect only logarithmic constants and disappear from letter classes. In particular, φ−1 and 1−φ give the same letter class despite their opposite signs.

Step 3 correctly substitutes the rational scalar coefficients as well as the tensor letters and keeps ζ(2), ζ(3) at length zero. It has one literal wording error: 𝔐 itself is not asserted by its displayed formula to be a sum with **only rational constant coefficients q**, since it includes the zeta multiples.

**Exact repair at line 262:** replace its first sentence by:

> The function R₃ is a sum of terms q c^e H_w(c), q∈ℚ, and 𝔐=R₃−ζ(3)Φ+ζ(2)Ψ, with Φ=1/2−(1−c²)(H₁+H₋₁)/(4c) and Ψ of length at most two.

The rest of Step 3 is unchanged. In particular, the zeta multiples cannot contribute length four.

### (iii) Naturality of δ₂,₂ — CORRECT

**Location:** `sec_g8_M3func.tex:267–271`.

On a pure tensor the identity is simply

δ₂,₂(fv₁⊗fv₂⊗fv₃⊗fv₄)=(fv₁∧fv₂)∧(fv₃∧fv₄).

This equals (Λ²Λ²f)δ₂,₂(v₁⊗v₂⊗v₃⊗v₄). Linear extension and substitution of coefficients give exactly line 269. Rational prefactors are scalars in this length filtration; no derivative of such a prefactor belongs in this operation.

### (iv) Injectivity, and preservation of the nonzero component — CORRECT

**Location:** `sec_g8_M3func.tex:273`.

The proof clears rational denominators in a relation of letter classes and forms G(X)=X^{n₀}(X−1)^{n₁}(X+1)^{n₋₁}. Unique factorisation of rational functions implies that [G(φ)]=0 means G(φ) is constant: there cannot be a nonconstant rational function having no finite zero or pole. A nonconstant rational φ is transcendental over ℂ, so substitution ℂ(X)→ℂ(t) is injective. Therefore G itself is constant, whose orders at 0,1,−1 force all three exponents to vanish. This argument also handles the case where a relevant preimage is infinity; a proof relying only on finding finite zeros separately would require more care.

The next implication is valid, not an assumption: an injective linear map f between these finite-dimensional ℚ-spaces has a left inverse L. Applying Λ² twice gives (Λ²Λ²L)(Λ²Λ²f)=identity. Thus the needed map is injective, and the independent vectors p∧r and q∧r remain independent. Extending scalars to ℂ(t) preserves this independence. Both coefficient functions 1−1/φ and −(1+1/φ) are nonzero for nonconstant φ. Hence their indicated linear combination cannot vanish.

### (v) Classical competitor and the final contradiction — CORRECT

**Locations:** `sec_g8_M3func.tex:177–210`, `:248`, `:275`.

The alphabet Σ_t includes the competitor's zeros/poles and coefficient poles through E(C), in addition to the pullback alphabet. The lemma is explicitly reapplied **in t**, over ℂ(t), on D′, with one common base point. Hyperlogarithmic independence applies on that disc even if C originally had factors of arbitrarily high weight. The classical-vanishing lemma concerns its degree-four component at that base point, so cancellation of higher lengths offers no loophole. The pullback has length at most four, making its degree-four component invariant under the base point. The contradiction is valid.

## 3. Exponentials, simultaneous sine/cosine, and conjugation — CORRECT

**Locations:** `sec_g8_M3func.tex:13–20`, `:239`, `:244`; omission acknowledged in `REPORT_FIX_G8BC.md:53`.

The stated class directly contains Li_n(t²), Li_n(t), and rational expressions in c=(t+t^{-1})/2 and s=(t−t^{-1})/(2i), since t=c+is. The elementary exponentials t and t² have no zeros, poles, or value 1 on the arc in question; branches exist on sufficiently small neighbourhoods. Other rational arguments with finitely many exceptional points can be handled by restricting to a smaller angle interval.

Conjugation and Re/Im are not listed as constructors in the definition. Nevertheless, **finite expressions obtained by adding these operations are already locally covered on the unit circle**, because they can be rewritten in the stated class. This needs a short explanation, not a stronger nonclassicality theorem.

Here is an exact one-line extension suitable after line 244:

> This also excludes finite expressions using conjugates, real parts or imaginary parts of such classical expressions along the unit circle: locally, if F is classical then F♯(t):=overline{F(1/overline{t})} is classical with rational arguments g♯(t)=ḡ(1/t) and coefficients h♯(t)=h̄(1/t), while on |t|=1 one has Re F=(F+F♯)/2 and Im F=(F−F♯)/(2i).

Here the bars on g,h conjugate their coefficients. For a branch, F♯ means its **reflected branch**, not an arbitrary choice of principal branch. Reflection preserves the defining exponential relation for P₁ and every differential relation P_k′=P_{k−1}g′/g, including the chain-rule factor from 1/t. Restrict to a small disc contained in both the original and reflected domains. Thus the extension applies to arbitrary allowed branch chains as well as usual polylogarithm branches.

For the specific principal-branch Clausen example, on 0<α<π,

Cl₂(α)=[Li₂(t²)−Li₂(t^{-2})]/(2i).

This is a restriction of a classical expression in t. Consequently expressions involving Cl₂ in the stated finite algebraic combinations are excluded as candidates for M₃. This does not assert that every real-analytic operation on a classical expression belongs to the class.

## 4. Independent exact pullback calculation — CORRECT

**Location of the checked claim:** `sec_g8_M3func.tex:280`; also `REPORT_FIX_G8BC.md:13`.

My code independently obtains

c=(t−i)(t+i)/(2t),

1−c=−(t−1)²/(2t),

1+c=(t+1)²/(2t).

For rows (0,i,−i,1,−1) and columns (x,y,z), the letter matrix is

```text
[-1 -1 -1]
[ 1  0  0]
[ 1  0  0]
[ 0  2  0]
[ 0  0  2]
```

The minor on rows (i,1,−1) has determinant 4. The ranks of φ_*, Λ²φ_*, and Λ²Λ²φ_* are all **3**. I computed the exterior maps using determinant minors, independently of the supplied scripts' implementation.

I also parsed the 56 weight-four terms directly from `sec_g8_M3alpha.tex:11–29` and verified equality with the supplied JSON terms. Reversal and H₁ signs yield 24 top tensors. Applying the determinant form of δ₂,₂ directly to their pulled-back tensor factors, with coefficients c→φ(t), agrees exactly with

(1−1/φ)(Λ²Λ²φ_*)(p∧r)−(1+1/φ)(Λ²Λ²φ_*)(q∧r).

There are **21 nonzero rational-function components out of 45**. For example, in the ordered basis above, the component along (e₀∧e_i)∧(e₀∧e₁) is −8t/(t²+1), already a nonzero witness. The output file lists all 21 components. The count is a count in this basis over ℂ(t), not a statement about specialising the letters to numbers at a fixed angle.

## 5. Title, “In words”, and “What is not claimed” — CORRECT AFTER REPAIR

**Locations:** `sec_g8_M3func.tex:1–3`, `:38–41`, `:244`, `:317`; `sec_g8_M3alpha.tex:234`.

The title and “In words” accurately describe the proved rational-in-t class. The interpretation of “needs depth-two functions” is expressly limited to failure of depth-one expressions; it does not prove depth two suffices or that the two-term remainder is classical. The disclaimers about a single angle, Goncharov's converse, and a classical remainder are appropriate.

The algebraic/genus boundary is inaccurate. Line 41 gives Li_n(√c) as an example of something outside the result and associates it with positive genus. But **√c alone is covered** by the general corollary: put c=u² and use Γ in the positive u-axis. The field ℂ(c,√c)=ℂ(u) has genus zero. Another counterexample to the alleged boundary is t=e^{iα/4}, explicitly covered in line 244 but not rational in e^{iα/2}.

The different assertion about the **joint** field ℂ(e^{iα/2},√c) is correct. Writing v=√c gives t²−2v²t+1=0 and (t−v²)²=v⁴−1. The squarefree quartic has four branch points, so its smooth projective curve has genus one. This does not make √c alone a positive-genus argument.

**Exact replacement for item (2) at line 41:**

> (2) The proof applies whenever c and all coefficient and argument functions in a proposed expression are rational functions of one common parameter t, with c nonconstant in t. Thus it covers √c alone via c=t² as well as the angle field ℂ(e^{iα/2})=ℂ(cos(α/2),sin(α/2)). It makes no general claim for algebraic arguments and coefficients whose joint function field with c has positive genus. For example, it does not treat arbitrary simultaneous dependence on e^{iα/2} and √c; their joint field has genus one. Panzer's statement [the existing citation] allows algebraic prefactors; the proof here uses prefactors rational in the chosen parameter.

**Exact replacement at line 317:** replace

> algebraic arguments beyond those rational in e^{iα/2}, that is, of positive genus

by

> general algebraic arguments and coefficients without a common rational parametrisation together with c

**Exact replacement in `sec_g8_M3alpha.tex:234`:** replace

> arguments algebraic over ℚ(c) beyond that (positive genus)

by

> general algebraic arguments and coefficients whose joint function field with c has positive genus

These changes correct an underdescription and misclassification of scope; they do not weaken the corollary or furnish a classical representation of M₃.

## 6. Revised W₃ definition and counts — CORRECT

**Locations:** `sec_g8_W3int.tex:14–37`, `:44–53`, `:134–153`, `:199–208`, `:217`, `:234`; corresponding summary in `sec_g8_M2alpha.tex:15`, `:26–35`.

The regularised recursive definition is consistent. I(0^r;y)=log^r(y)/r! fixes the tangential normalisation. Every other word contains a nonzero letter. At the first nonzero letter outside a pure-zero tail one integrates log^r(y) against a bounded kernel, obtaining O(y(1+|log y|^r)). Subsequent zero letters divide by y but preserve local integrability after integration; subsequent nonzero letters improve the power of y. This proves the needed remaining improper integrals converge absolutely. The nonzero letters are negative and the endpoint is positive, so there is no further singularity along the path.

For a,b not both zero, both regularised two-letter terms vanish at 0 and their sum has derivative d(λ_aλ_b); hence their sum is λ_aλ_b. For a=b=0 the same identity is the prescribed normalisation 2I(0,0)=log²y. Integrating once more proves the claimed length-three expansion. Thus no hidden extra constants are introduced by the shuffle step.

Independent coefficient extraction from the quadratic forms in λ gives:

| Quantity | Result |
|---|---:|
| Nonzero integrals of products λ_aλ_b, a≤b | 64 |
| Counts for outer index j=0,1,2,3,4 | 10, 11, 13, 15, 15 |
| Diagonal / off-diagonal products | 23 / 41 |
| Distinct ordered words after collecting diagonals | 105 |
| Uncollected summands | 128 |
| Words ending in 0 | 21 |
| Words (j,0,0) | j=1,2,3,4 |
| Pure word (0,0,0) | absent |

Multiplication of every coefficient by c gives a polynomial over ℚ of degree at most two. This confirms ℚ+ℚc+ℚ/c for both product coefficients and ordered-word coefficients. These are formal expansion counts, not dimensions or minimality claims. They remain correct despite potential identities between integrals or special parameter values.

As an additional implementation check, I verified all 64 shuffle-plus-integration identities in exact power/log series through y⁷ at c=3/5. This finite series check supplements, and does not replace, the preceding all-orders derivative/boundary argument.

I also recomputed all three partial-fraction identities of N1, the subtraction F_A−V′, and V(dlog g/dt)=1/[2(t+1)] as exact rational identities. Independent 40-digit quadratures of the differentiated defining integrand and W₃♯ give |M₂−W₂−W₃♯|=2.30×10⁻⁴¹ at c=3/5 and 5.74×10⁻⁴² at c=√3/2. In the latter case W₃♯=−0.376104091935813048896981…, agreeing with the supplied table. Numerical agreement is a check, not the proof.

## 7. New errors or inconsistencies in the four sections — CORRECT AFTER REPAIR

**Locations:** those listed in items 2 and 5 above; `sec_g8_M2alpha.tex:26–35`; `sec_g8_M3alpha.tex:44`, `:62`, `:126`, `:141`, `:234`; `sec_g8_W3int.tex:32–53`, `:206`.

The defects identified in the changed/new text are:

1. Step 1 changes “coefficient poles” to “poles of C”; use the original finite excluded set.
2. Step 3 describes all of 𝔐 as having rational constant coefficients despite its zeta multiples; distinguish R₃ or allow complex constant coefficients.
3. The algebraic/genus disclaimers misclassify √c and fail to reflect the general rational-parametrisation corollary; repair all three occurrences specified in item 5.

I found no further mathematical error in the revised W₃ material or the M₂ summary of it. The revised M₃ cross-references correctly point to the function-level invariant and angle-field corollary, apart from the genus disclaimer. The newly added conjugation explanation in item 3 is recommended for clarity; its absence does not invalidate the stated result.

There is no pre-fix snapshot in this review folder, so I cannot certify which other lines changed during the fix pass. “New” here refers to the corollary and the revised areas identified in the fixer's report, not to a verified historical diff. Several external sections and cited check scripts/outputs are also absent from this extracted folder; that limits build and historical reproducibility, not the exact calculations performed above.

## One sentence the paper may truthfully say

On no nonempty open interval J⊂(0,π) can M₃(α) be represented by a finite ℂ(e^{iα/2})-linear combination of products of complex constants, branches of logarithms and classical polylogarithms of arbitrary weights with arguments rational in e^{iα/2} (equivalently rational in cos(α/2) and sin(α/2)), including finite conjugate and real/imaginary-part combinations interpreted locally as above.

## What I did not check

- I accepted the previously refereed identification of M₃ with R₃−ζ(3)Φ+ζ(2)Ψ and the earlier all-n nested-sum proof. I did independently verify the supplied 56-term TeX/JSON match, its pullback and its obstruction; I did not rederive the 74-term expansion from the original corner integral.
- I did not recheck the full two-variable Li_{2,2} calculations, the constant-letter remainder, or the “one term cannot” theorem beyond reading their interface with this fix.
- I did not rerun the full M₂ classical-form/all-n machinery or prove its displayed HPL formula anew. My numerical W₃ check used the differentiated coefficient integrand, not a Cauchy extraction from M(u;α).
- I did not validate the historical timings, stored outputs, earlier referees' work, peer-review status, or the claimed complete-wrapper build. This folder lacks some wrapper inputs and several cited scripts/outputs.
- I make no claim about individual numerical angles, general positive-genus algebraic arguments/coefficient fields, sufficiency of depth two, a classical remainder, or the motivic identification of the explicit tensor map.

The check's final hash comparison confirms the supplied source/data files remained unchanged during execution; the environment's live `run.log` is excluded from that comparison.
