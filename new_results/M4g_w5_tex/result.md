# M4g — LaTeX for the W5 / M4 result: supplement, checks, paper fragment

## File list (in `new_results/M4g_w5_tex/`)

| file | role |
|---|---|
| `gen_supplement.py` | generic generator: `gen_supplement.py [TERMS_JSON [OUT_TEX]]`; builds the table **only** from the JSON |
| `supplement_W5.tex` / `.pdf` | the supplement (from `R5_nice_reduced.json`), standalone `article` + amsmath/amssymb/geometry |
| `supplement_W5_best.tex` / `.pdf` | same generator run on M4f's `R5_best.json` (70 terms, incl. `Li5(b⁴)`, `Li5(1−b⁴)`) |
| `supplement_W5_old.tex` | saved pre-generic copy (for the diff below) |
| `check_supplement.py` | round-trip parser + 400-digit evaluator: `check_supplement.py [TERMS_JSON [TEX_PATH]]` |
| `paper_fragment.tex` | 9 source lines + one `multline*` display; `\input`-ready (planner-edited wording) |
| `fragment_test.tex` | minimal `article` + amsmath wrapper that `\input`s the fragment |
| `progress.log` | step log |

## Round-trip and 400-digit check (printed output, `.venv` python)

On the default input (`R5_nice_reduced.json`):
```
round-trip (a): all 68 classical terms match the JSON exactly
round-trip (a): all 4 N5 terms match cobracket.json exactly
classical part evaluated: 37.653188
W5 estimate : -0.18235662466773630013
W5 (5400bit): -0.18235662466773630013
difference  : -9.0815e-427
CHECK (b): |W5est - W5| < 1e-390  PASS (got ~9.08e-427)
```
On `../M4f_w5_extras/R5_best.json` (exists; was run):
```
round-trip (a): all 70 classical terms match the JSON exactly
round-trip (a): all 4 N5 terms match cobracket.json exactly
difference  : 2.982e-427
CHECK (b): |W5est - W5| < 1e-390  PASS (got ~2.98e-427)
```
Method (b): direct series for `Li_k` and `Li_{3,2}` (running inner sum), `Li5(−1)=−(15/16)ζ(5)`, `L(χ8,5)=8^{−5}[ζ(5,1/8)+ζ(5,7/8)−ζ(5,3/8)−ζ(5,5/8)]` via Hurwitz ζ, `mp.mp.dps=430`. (The earlier spurious `~3e-12` residuals were a bug: `mp.dps=430` sets a module attribute and leaves 15-digit precision; fixed to `mp.mp.dps`.)

## Diff of the regenerated supplement vs the accepted copy (only the three removed leading `+`; Li4/w3w2 first lines start with `−` and are unchanged)
```diff
< N_5={}&+\bigl(779-\tfrac{695}{56}\sqrt2\bigr)\,\operatorname{Li}_{3,2}\bigl(b,1\bigr)\\
> N_5={}&\bigl(779-\tfrac{695}{56}\sqrt2\bigr)\,\operatorname{Li}_{3,2}\bigl(b,1\bigr)\\
< W_5-N_5={}&+\bigl(\tfrac{391577}{252}+\tfrac{1001239}{2016}\sqrt2\bigr)\,\operatorname{Li}_{5}\bigl(b\bigr)\\
> W_5-N_5={}&\bigl(\tfrac{391577}{252}+\tfrac{1001239}{2016}\sqrt2\bigr)\,\operatorname{Li}_{5}\bigl(b\bigr)\\
< &+\bigl(\tfrac{1692679}{92160}+\tfrac{1026973}{2580480}\sqrt2\bigr)\,\lambda^{5}\\
> &\bigl(\tfrac{1692679}{92160}+\tfrac{1026973}{2580480}\sqrt2\bigr)\,\lambda^{5}\\
```

## Compile results (`/opt/homebrew/bin/tectonic -X compile <file>`)

| file | tectonic tail line | errors | overfull hbox >10pt | pages |
|---|---|---|---|---|
| `supplement_W5.tex` | `note: Writing `​`supplement_W5.pdf`​` (41.66 KiB)` | 0 | 0 | **3** |
| `supplement_W5_best.tex` | `note: Writing `​`supplement_W5_best.pdf`​` (41.78 KiB)` | 0 | 0 | **3** |
| `fragment_test.tex` | `note: Writing `​`fragment_test.pdf`​` (23.49 KiB)` | 0 | 0 | **1** |

(Each run: "Running TeX … Running xdvipdfmx … Writing …pdf"; no `warning:` lines. An earlier 274.9pt overfull from the 50-digit `W_5`/`M_4` strings was fixed by moving them into a `gather*`.)

## The fragment, verbatim (`paper_fragment.tex`, planner-edited)
```tex
We come to $M_4$. Splitting by weight, $M_4=W_5+W_4''+W_3^{(4)}+T$, and the last three are already in closed form (proved; $W_4''$ by an exact computer-assisted certificate, M4b2): $W_4''$ in tetralogarithms, $W_3^{(4)}=W_3'-\frac18\zeta(3)-\frac1{24}\pi^2\lambda-\frac1{12}\lambda^3$, and $T=(\pi^2+6\lambda^2)/24$. What remains is $W_5$, the part of weight five.
$W_5$ is not classical at the level of symbols; for the number itself this assumes the standard period conjecture. With $X=\lambda$ and $Y=L$, the antisymmetrised $(3,2)$ part of its symbol pairs to $84$ against $[X,[X,Y]]\otimes[X,Y]$ and to $-38\sqrt2$ against $[Y,[X,Y]]\otimes[X,Y]$. Every $\operatorname{Li}_5$ and every product of lower-weight constants pairs to zero, so no combination of $\operatorname{Li}_5$'s and products has the symbol of $W_5$.
The $(3,2)$ cobracket fixes a depth-two correction exactly. With $\operatorname{Li}_{s,t}(x,y)=\sum_{n>m>0}x^ny^m/(n^sm^t)$,
\begin{multline*}
N_5=\Bigl(779-\frac{695\sqrt2}{56}\Bigr)\operatorname{Li}_{3,2}(b,1)-672\,\operatorname{Li}_{3,2}\Bigl(\frac1{\sqrt2},1\Bigr)\\
-(374-76\sqrt2)\operatorname{Li}_{3,2}(-b,1)-\Bigl(405-\frac{4951\sqrt2}{56}\Bigr)\operatorname{Li}_{3,2}(b,-1).
\end{multline*}
The difference $W_5-N_5$ is a combination, with coefficients $p+q\sqrt2$ ($p,q$ rational), of twelve $\operatorname{Li}_5$'s, $\zeta(5)$, $L(\chi_8,5)$ and products of lower-weight constants. The twelve arguments are $b$, $-b$, $-b^2$, $1/\sqrt2$, $\frac12$, $1-1/\sqrt2$, $2-\sqrt2$, $2b$, $(2-\sqrt2)/4$, $(2+\sqrt2)/4$, $-b/2$ and $(4-3\sqrt2)/8$.
The identity holds numerically to 1596 digits, with exact agreement of the symbol and the $(3,2)$ cobracket. The coefficient of $L(\chi_8,5)$ depends on the representative, because $L(\chi_8,5)$ lies in the $\operatorname{Li}_5$ span. The full table is in the Supplement. $M_4$ is the first squircle coefficient that needs a depth-two polylogarithm (at the level of symbols).
```

## Rerun commands (from `new_results/M4g_w5_tex/`)
```sh
PY='<project>/.venv/bin/python'
nice $PY gen_supplement.py                                   # -> supplement_W5.tex   (default: R5_nice_reduced.json)
nice $PY check_supplement.py                                 # round-trip + 400-digit check
nice $PY gen_supplement.py ../M4f_w5_extras/R5_best.json supplement_W5_best.tex
nice $PY check_supplement.py ../M4f_w5_extras/R5_best.json supplement_W5_best.tex
nice /opt/homebrew/bin/tectonic -X compile supplement_W5.tex
nice /opt/homebrew/bin/tectonic -X compile supplement_W5_best.tex
nice /opt/homebrew/bin/tectonic -X compile fragment_test.tex
```

## Notes
- The generator is fully data-driven: group orders are built from the keys present (W4's seven `Li5` arguments first, then `b², −1/√2, 2b, (2−√2)/4, (2+√2)/4, −b/2, (4−3√2)/8, b⁴, 1−b⁴`, then `ζ(5)`, `L(χ8,5)`; Li4/w3w2/pure templates filtered to present keys; every key asserted placed exactly once). All 18 argument forms are checked numerically against the JSON field elements at generation time.
- The supplement states: numerical identity only (1596 digits vs certified W5, exact symbol and (3,2) cobracket), not a proof; the `L(χ8,5)` coefficient is representative-dependent because `L(χ8,5)` lies in the `Li5` span; `W_5` and `M_4` are displayed to 50 digits each.