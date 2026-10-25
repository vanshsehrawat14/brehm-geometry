# Reference specification

What `brehm_cpp` computes, what a certificate asserts, how each checker establishes it, and what is
not established. Written from Brehm's paper (U. Brehm, *Extensions of distance reducing mappings to
piecewise congruent mappings on R^m*, J. Geom. 16 (1981) 187–193, doi:10.1007/BF01917587), not
from the code.

## 1. Supported input

| Item | Supported | Rejected (exit 3) |
|---|---|---|
| paper `D` | bounded convex polygon with positive area, rational vertices in boundary order (either orientation, collinear vertices allowed) | non-convex, self-intersecting, zero area |
| pairs `(a_i, b_i)` | `n ≥ 1` rational points of `R²`, `a_i ∈ D` (boundary allowed), `|b_i − b_j| ≤ |a_i − a_j|` exactly; repeated sources with equal targets | expanding pairs (including one source with two targets), sources outside `D` |
| order | any permutation (default: as listed) | not a permutation |
| `g` | linear part of Brehm's initial isometry, any exactly orthogonal rational 2×2 matrix (rotations and reflections with rational entries) | not orthogonal |
| `L` | rational, `L² ≥ 2 max_i |a_i − a_{σ(1)}|²` (Brehm's `d(x_i, x_1) ≤ r` with `r = L/√2`) and `D ⊂ a_{σ(1)} + [−3L, 3L]²`; default `2·max(1, L¹ radius of D about a_{σ(1)})` | below either bound |
| numbers | JSON numbers, or strings holding an integer, `p/q` (`q > 0`) or a decimal with exponent `|e| ≤ 1000`, read exactly (§3) | anything else, including spaces, `1/0`, `NaN` |

Not supported: dimensions other than 2 → 2, irrational data (for example rotations by angles with
irrational cosine), unbounded or non-convex paper. Nothing outside this table is claimed.

## 2. The map being computed

For the pairs in insertion order `(a_1, b_1), …, (a_n, b_n)`, Brehm's construction with these
choices defines maps `f_1, …, f_n : R² → R²`:

* `f_1(x) = b_1 + g·φ(x − a_1)`, `φ` acting on each coordinate by `t ↦ t` on `[−L, L]`,
  `2L − t` for `t ≥ L`, `−2L − t` for `t ≤ −L` (Brehm p. 188, with `r√m = L`).
* For `i ≥ 1`: `U_i = {x : |x − a_{i+1}| < |f_i(x) − b_{i+1}|}`; `f_{i+1} = f_i` off `U_i`;
  on `U_i`, `f_{i+1}` is the unique distance-reducing extension (Brehm, Remark 2, p. 192), realized
  by the cones `conv((K ∩ ∂U_i) ∪ {a_{i+1}})` over the cells `K` of `f_i`.

`B(D, data, order, g, L) := f_n|_D` is the object `brehm_cpp construct` claims to output. It depends
on the order, on `g` and on `L` (see §7).

## 3. Certificate (`brehm-rational-v1`)

JSON with `domain`, `sources`, `targets` (exact numbers as strings) and `pieces`: a list of cells
`{polygon, A, b}` meaning `x ↦ A x + b` on the closed convex polygon (`name` fields are optional).
`brehm_cpp` adds a `construction` record `{order, g, L}` that the validity checkers ignore and the
evaluator (§6) uses.

**Portable reading rules** (`json.hpp`, `ref_check.py`, `brehm_eval.py`; the figure
tool `brehm_figure.py` reads with `ref_check.py`'s reader):

* a number is a JSON number or a string matching `[+-]?(D/D | (D.D* | .D)([eE][+-]?D)?)` with `D` a
  run of decimal digits, no spaces, nonzero denominator, and explicit exponent `|e| ≤ 1000` (so a
  short literal cannot demand a gigabyte-sized number);
* duplicate keys in an object and nesting deeper than 64 levels are rejected. Python's `json` keeps
  the last duplicate and a naive reader the first, so without this rule two checkers could accept two
  different certificates from one file;
* use standard JSON numeric tokens and ASCII digits in rational strings. The C++ and Python
  readers have acceptance differences outside this portable subset, including Unicode escapes
  and nonstandard numeric tokens; see [LIMITATIONS.md](LIMITATIONS.md).

## 4. Validity: conditions R1–R6

A certificate is **valid** if

* **R1** `D` is a convex polygon with positive area;
* **R2** every cell `K_i` is a convex polygon with positive area;
* **R3** `K_1 ∪ … ∪ K_m = D`;
* **R4** every `A_i` is orthogonal: `A_iᵀA_i = I` exactly;
* **R5** `A_i x + b_i = A_j x + b_j` for every `x ∈ K_i ∩ K_j`, for all `i, j`;
* **R6** every `a_k ∈ D`, and `A_i a_k + b_i = b_k` for every cell containing `a_k`.

**Theorem V.** If R1–R6 hold, `F(x) := A_i x + b_i (x ∈ K_i)` is a well-defined, continuous map
`D → R²`, an affine isometry on each `K_i`, with `F(a_k) = b_k`; it is piecewise congruent in Brehm's
sense and distance-reducing.

*Proof.* R5 makes `F` well defined; it is continuous by the pasting lemma (finitely many closed sets
covering `D`, R3). R4 makes each piece an isometry, R6 gives `F(a_k) = b_k`. For `x, y ∈ D` the
segment `[x, y]` lies in `D` (R1, convex) and meets finitely many convex cells in finitely many
subsegments; `F` preserves the length of each, and by continuity `F([x, y])` is a path of length
`|x − y|`, so `|F(x) − F(y)| ≤ |x − y|` (Brehm, Remark 1, p. 188). ∎

R1–R6 is the minimal definition: it does not require disjoint interiors or a conforming complex.
All three checkers verify something at least as strong.

## 5. The three checkers, and an audit of the earlier specification

| Checker | Coverage (R3) | Continuity (R5) | Code |
|---|---|---|---|
| repository | cells in `D`, pairwise interiors disjoint, `Σ area = area(D)` | maps agree at both endpoints of every shared-edge overlap; vertex-only contacts by the sector walk | `python/core/exact_certify.py` |
| C++ | the 1-chain `Σ ∂K_i − ∂D` vanishes on every elementary interval, with exactly two edges on each | maps agree at both endpoints of every interior elementary interval; vertex contacts by the cyclic order of cells | `check.cpp` |
| reference | vertical slabs: in each open slab between consecutive vertex x-coordinates, cells stack from `D`'s lower line to its upper line with equal lines at every contact | for every pair of cells with intersecting bounding boxes, maps agree at all candidate extreme points of `K_i ∩ K_j` | `tools/ref_check.py` |

**Soundness of each coverage test.**

* *Area* (repository): a point of `D` outside `∪K_i` has a neighbourhood in `D` of positive area
  missed by the closed cells; with disjoint interiors `area(∪K_i) = Σ area = area(D)`, a
  contradiction. Containment of each cell in `D` is checked directly.
* *Boundary cancellation* (C++): `w(x) = Σ_i wind(∂K_i, x) − wind(∂D, x)` is the winding number
  of the 1-cycle `Σ ∂K_i − ∂D`; if that cycle is zero as a chain, `w ≡ 0` off the edges. Each
  counter-clockwise convex cell has winding number 1 inside, so a.e. point of `D` lies in exactly one
  cell and a.e. point outside in none; closedness gives R3 and `K_i ⊂ D`.
* *Slabs* (reference): no vertex has its x-coordinate inside an open slab, so each cell meeting the
  slab spans it and is bounded there by one lower and one upper line; equal lines at every contact
  make the stack tile `D ∩ slab` for every x in the slab, not only at the midline. The slabs cover `D`
  up to finitely many vertical lines; closedness gives R3.

**Soundness of each continuity test.** Two affine maps agreeing at the extreme points of a convex set
agree on it. The reference checker tests every pairwise intersection directly (extreme points of
`K_i ∩ K_j` are vertices of one polygon in the other or edge crossings), so it needs no argument
about vertex-only contacts. The other two test shared edges only and close the gap with the cyclic
(or, on `∂D`, linear) order of cells around a point, which is valid because the cells tile a
neighbourhood of the point in `D`.

**Audit of the earlier specification** (`python/docs/certificate_lemma.md` of the development
repository, not included here), which both earlier checkers were written from. Its hypotheses 1–4 imply R1–R6: hypothesis 1 gives R2–R3 by the area argument,
2 = R4, 3 together with the sector walk gives R5, 4 = R6 (sources in the union of pieces, which is
`D`). The proof is correct. It is stronger than necessary (it needs disjoint interiors) and states the
vertex-contact step only in prose. Nothing in it is wrong for rational certificates. Its remark on
floating-point certificates does not apply to `brehm_cpp`, which is exact.

**Independence.** The reference checker was written from R1–R6 above, uses a third coverage method
and a direct pairwise continuity test, and imports nothing from the project. It shares the author
lineage and the JSON format with the other two. An independent human reader of §4–§5 is still needed.

**Shared assumptions (audit, 16 October).** What a common error could slip through:

1. *The specification.* All three implement R1–R6, so an error in §4 would be shared. Theorem V is
   short and uses only the pasting lemma and the convexity of `D`.
2. *The reader.* Until 16 October the C++ reader kept the first of two duplicate keys and Python the
   last, so one file could present two different certificates to different checkers; numbers such as
   `1 / 2`, `1_0` or `1/-2` were read differently. All readers reject duplicates; remaining grammar differences are documented in §3 and
   LIMITATIONS.md, and `tests/malformed.py` runs every tool on the same hostile files. Until the
   release checks of 16 October the 64-level nesting limit was enforced only in C++: the Python tools
   accepted a file nested 200 levels deep in an ignored key, and on Python 3.11 the evaluator crashed
   (stack overflow) on 100,000 levels. Both Python tools now count levels the way `json.hpp` does,
   before and after parsing, and the tests probe 64 and 65 levels. The
   repository checker (`python/core`) has no reader of its own here: `run_tests.py` feeds it what
   `ref_check.py` read.
3. *Arithmetic.* All exact (GMP `mpq_class`, Python `Fraction`); no tolerance anywhere.
4. *Data binding.* A checker verifies the certificate's own domain and pairs. Only
   `ref_check.py --instance` (domain, pairs) and `brehm_eval.py --instance` (pairs, order, `g`, `L`)
   tie them to the instance; `brehm_cpp check` does not.
5. *Lineage.* One human author directed development with Claude Code, OpenAI Codex and ChatGPT.
   These are complementary implementations, not independent human reviews; see README's disclosure.
6. *Availability.* The repository checker is `python/core/exact_certify.py`, shipped here. Without
   `python/core`, `run_tests.py` reports two checkers.

## 6. Faithfulness to Brehm, without cells

**The transcription.** §2 was checked against the paper (pp. 188–192) on 16 October:
`f_1 = g ∘ h(x − x_1)` with `h` the coordinatewise `φ` (p. 188); `C = x_1 + [−3r√m, 3r√m]^m`
(p. 189); `U` as in §2 (p. 190); `U ⊂ C°` (p. 191, item 4); cones
`K̃_ν = conv((K_ν ∩ ∂U) ∪ {x_{i+1}})` (p. 192). The paper has two remarks numbered 2: on p. 188 (a
piecewise congruent map on a compact convex polytope is the only distance-reducing map with its
values at the cell vertices) and on p. 192 (`f_{i+1}` is the unique distance-reducing extension of
`f̃_i`).

**What is taken from Brehm.** For short prescribed data, exactly orthogonal `g`, and admissible `L` (`d(a_i, a_1) ≤ r = L/√2`):

* (B1) every `f_i` is piecewise congruent (his Theorem), hence distance-reducing on `R²` (Remark 1,
  p. 188);
* (B2) `f_{i+1} = f_i` off `U_i`, and `f_{i+1}(a_{i+1}) = b_{i+1}` (p. 192);
* (B3) `U_i ⊂ C°` (p. 191).

**Lemma E1.** Let `f = f_i`, `(a, b) = (a_{i+1}, b_{i+1})`, `x ∈ U_i`, and `p = a + t(x − a)` with
`t > 1` and `|f(p) − b| = |p − a|`. Then `f_{i+1}(x) = b + (f(p) − b)/t`.

*Proof.* `p ∉ U_i`, so `f_{i+1}(p) = f(p)` (B2); `f_{i+1}(a) = b`. `f_{i+1}` is 1-Lipschitz (B1) and
`|f_{i+1}(p) − f_{i+1}(a)| = |p − a|`, so for `x` on `[a, p]` both inequalities
`|f_{i+1}(x) − b| ≤ |x − a|` and `|f_{i+1}(x) − f(p)| ≤ |p − x|` are equalities, which forces
`f_{i+1}(x)` to be the point of `[b, f(p)]` at distance `|x − a|` from `b`. ∎

Such a `p` exists: `t ↦ |f(a + t(x − a)) − b| − t|x − a|` is nonincreasing (B1), positive at
`t = 1`, and `≤ 0` where the ray leaves `C` (B3). On each affine piece of `f` along the ray the
squared discrepancy is linear in `t`, so the exit parameter is rational.

**Lemma E3 (the evaluator returns only correct values).** Every value `e_k(x)` returned by
`tools/brehm_eval.py` equals `f_k(x)`.

*Proof.* Induction on `k`. `e_1` is the formula for `f_1`. For `k > 1` the evaluator computes
`e_{k−1}(x)`, correct by induction, and decides `x ∈ U_{k−1}` by the exact test
`|x − a_k|² < |e_{k−1}(x) − b_k|²`. If `x ∉ U` it returns `e_{k−1}(x)` (B2); if `x = a_k` it returns
`b_k` (B2). Otherwise it returns `b_k + (e_{k−1}(p) − b_k)/t` only after checking exactly that
`|e_{k−1}(p) − b_k|² = t²|x − a_k|²` and establishing `t > 1` (directly
for search proposals, or from (B3) at the protected-square boundary), with `e_{k−1}(p)` correct by induction: that is E1.
The search that proposes `p` (affine germs, a bracket that ends on `∂C`, bisection) decides only
which `p` is tried and how long that takes. ∎

So a wrong germ or any other search heuristic can make the evaluator slow or undecided (exit 4), never
wrong for admissible data and choices. Because (B1)–(B3) need admissibility, the evaluator checks
that `g` is exactly 2×2 with `gᵀg = I`, using rational arithmetic, and checks
`L² ≥ 2 max |a_i − a_1|²` itself. It rejects inadmissible choices with exit 3. Both rotations and
reflections are allowed. Validity and instance-domain binding still require a validity checker.

**Theorem E2.** Let a certificate satisfy R2–R4 and let `T_K(v) = f_n(v)` for every cell `K` and
every vertex `v` of `K`, where `T_K` is `K`'s own isometry. Then `F = f_n` on `D`.

*Proof.* Let `x = Σλ_j v_j ∈ K` (convex weights; `K` is the hull of its vertices by R2). For any point
`y` and the points `w_j = T_K v_j`, whose `λ`-barycenter is `T_K x`,
`Σλ_j |y − w_j|² = |y − T_K x|² + Σλ_j |T_K x − w_j|²`. Take `y = f_n(x)`: by (B1),
`|f_n(x) − w_j| = |f_n(x) − f_n(v_j)| ≤ |x − v_j| = |T_K x − w_j|`, so `|f_n(x) − T_K x|² ≤ 0`.
The cells cover `D` (R3). ∎ This is the argument of Brehm's Remark 2 on p. 188.

E2 does not assume that `f_n` is affine on `K`; it proves it. In particular:

* **Hidden creases.** A crease of `f_n` through the interior of a cell would make `f_n` non-affine on
  that cell, which E2 excludes once the vertices agree. Agreement at the vertices is therefore
  enough, and validity alone is not: `tests/run_tests.py` folds the tip of a corner cell over a short
  chord. The result is a valid extension (all three checkers accept it) that is not Brehm's map, and
  the evaluator rejects it at the tip's corner. Same test on four certificates.
* **Refinement boundaries.** Cells finer than Brehm's (a cell cut in two with one map) do not matter:
  E2 is per cell, and the evaluator accepts such a refinement (tested).
* **Degeneracies.** Collinear vertices are allowed (E2 needs only `K = conv(vertices)`). Vertices on
  `∂U_i`, on the fold lines of `f_1` and at data points, and rays whose zero-discrepancy set contains a segment, are
  decided by the exact tests of E3. A zero-area cell is excluded by R2 and would be harmless (the
  argument holds on any convex hull).

**What this proves and what it assumes.** "R1–R6 pass and the evaluator agrees at every vertex"
proves that the certificate equals Brehm's `f_n` for the recorded `(order, g, L)`, using neither
constructor. With `--instance`, `brehm_eval.py` also checks that the certificate's pairs and recorded
choices are the instance's (defaults resolved as `brehm_cpp` resolves them). Assumed: (B1)–(B3) from
the paper; the transcription in §2, checked against the paper but not yet independently; and that
`brehm_eval.py` implements E3. It is unverified Python, tested by 31 examples, 2,350 fuzz instances
and the mutation and attack tests.

**What Brehm's uniqueness statement establishes.** Remark 2 (p. 192) says each step is forced:
given `f_i` (hence `U_i`) and the new pair, `f_{i+1}` is the only distance-reducing extension of
`f_i|_{R² ∖ U_i}` plus `a_{i+1} ↦ b_{i+1}`. By induction `f_n` is determined by `f_1` and the order,
that is, by `(g, r, order)`. It does **not** say the final extension is unique for the data, and it
does not make the output independent of `g`, `r` or the order. All three are explicit choices; `examples/options_small_L.json`,
`options_order.json` and `options_g_rotation.json` illustrate varying them on generic data.

## 7. Termination and failure

| Tool | Terminates | Failure modes |
|---|---|---|
| `brehm_cpp construct` | **yes, structurally**: `n − 1` insertions, each visits every current cell once with a bounded number of exact operations per vertex; no data-dependent loops. At most `9·2^{n−1}` cells (each cell gives at most a kept part and a cone) | memory and time grow with cells and bit length (the cell bound is exponential; a tight worst-case bound is not established here); no limit is enforced. The only internal exception, a degenerate cone, is unreachable: it needs `a` on the bisector `u·x = c`, but `c − u·a = |u|²/2 > 0` whenever `u ≠ 0` |
| `brehm_cpp check`, `ref_check.py` | **yes**: finite loops over cells, edges, slabs and cell pairs | memory and time only |
| `brehm_eval.py` | **not proved.** The exit search starts from the bracket `[1, t_C]` (`t_C` where the ray leaves `C`; the discrepancy is `≤ 0` there by (B3)), so there is no unbounded doubling phase. Within the bracket it has a step budget of `400 + 8·(bit length of the numbers involved)`. On a piecewise affine `f_k` it succeeds once both bracket ends lie on the piece that holds the exit, which takes `O(bits)` halvings; that is an argument, not a proof | budget exhausted: **exit 4, "undecided"**, never a verdict. Seen once before the 16 October fix (a search for a hint at 5,000-digit coordinates); none since, in 31 examples, 2,350 fuzz instances and the attack tests |

**Exit codes** (all four command-line tools): 0 ok, accepted or agrees; 1 rejected (a certificate
that fails R1–R6, or disagrees with Brehm's map); 2 usage; 3 invalid input (unreadable file, number
grammar, duplicate keys, nesting, inadmissible options, non-convex paper, expanding pairs); 4 internal
error, I/O error, or undecided (evaluator); 5 out of memory (`brehm_cpp`: a GMP allocation hook and
`std::bad_alloc`; Python: `MemoryError`). Only exit 0 is a positive verdict. On Linux the kernel's
OOM killer can end a process first (SIGKILL), which no program can report. Parsing and diagnostic exceptions to this contract are documented in LIMITATIONS.md.
The current validation record is in VALIDATION.md.

**Correctness of the code is not formally proved.** It is established per output: R1–R6 by three
checkers, faithfulness by E2 and E3. Test results (examples, fuzzing, sanitizers) are evidence about
the code, not proofs about it.

## 8. Not covered

* Physical (layered) flat foldability. A valid certificate is a mathematical extension, read as a
  crease pattern it satisfies Kawasaki's condition at every interior vertex, but this software does not
  decide whether it has a flat folded state.
* Whether the input file is the data you meant. `ref_check.py --instance` compares the certificate's
  domain and pairs with an instance file.
* The insertion-order experiments `brehm_cpp orders-exhaustive` and `orders-anneal` (README). They
  report cell counts of `construct` for other orders, with `g = id` and the default `L` of the first
  inserted source; they write no certificate, and their counts are not checked.
