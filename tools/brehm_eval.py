"""Evaluate Brehm's extension pointwise, without cells, and test a certificate against it exactly.

Brehm (J. Geom. 16, 1981, pp. 188-192) sets f_1(x) = g(h(x - a_1)), h the coordinatewise fold at
+-L (his r*sqrt(m)), and for each later pair (a, b), with f = f_i,

    U = {x : |x - a| < |f(x) - b|},     f_{i+1} = f off U,

where f_{i+1} is the unique distance-reducing extension of f restricted to (R^2 minus U) plus a -> b
(his Remark 2, p. 192). This file never forms cones or cells. It uses one lemma (docs/SPEC.md, E1):

  Let x be in U and p = a + t (x - a) with t > 1 and |f(p) - b| = |p - a|. Then
  f_{i+1}(x) = b + (f(p) - b) / t.

  (p is not in U, so f_{i+1}(p) = f(p); f_{i+1}(a) = b; a 1-Lipschitz map with
  |f_{i+1}(p) - f_{i+1}(a)| = |p - a| maps the segment [a, p] linearly onto [b, f(p)].)

Such a p exists: t -> |f(a + t(x - a)) - b| - t|x - a| is nonincreasing, positive at t = 1, and
negative where the ray leaves Brehm's protected cube. The search for p uses the local affine piece
("germ") of f only to propose candidates; every candidate is accepted only after the exact test
|f(p) - b|^2 = |p - a|^2, so a wrong germ can slow the search down but never change a value.

Faithfulness (SPEC E2): if a certificate passes the R1-R6 check and its map agrees with f_n at every
vertex of every cell, it equals f_n on the whole domain (Brehm's Remark 2, p. 188: on each cell the
unique 1-Lipschitz extension of an isometric map on the vertices is the affine isometry).

    python brehm_eval.py CERT.json [--instance INSTANCE.json] [--L inf] [--points K --seed S]

The choices (order, g, L) are read from the certificate's "construction" record; without one, the
defaults of python/core/original_brehm.py are used (input order, g a translation, L = 2R with
R = max(1, L1 radius of the domain about the first source)). With --L inf the folds of f_1 are
dropped; the run then records the largest |x - a_1|_inf at which f_1 was queried (M), and the result
holds for every L > M (the evaluation never sees the folds). With --instance it also checks that
the certificate's pairs and its recorded choices are those the instance asks for (defaults
resolved as brehm_cpp resolves them). Standard library only.

Exit codes: 0 agrees, 1 disagrees, 2 bad usage, 3 unreadable input (including duplicate JSON keys
and numbers outside the grammar of SPEC §3), 4 undecided (no exit point found within the step
budget) or internal error, 5 out of memory. Only exit 0 is a positive verdict.
"""
import json
import random
import re
import sys
from fractions import Fraction

sys.setrecursionlimit(100000)

NUMBER = re.compile(r"[+-]?(?:\d+/\d+|(?:\d+\.?\d*|\.\d+)(?:[eE]([+-]?\d+))?)")


def num(v):
    """An exact number: a Fraction, or a string that is an integer, a decimal with |exponent| <= 1000,
    or p/q with q > 0 (the grammar of brehm_cpp and tools/ref_check.py)."""
    if isinstance(v, Fraction):
        return v
    m = NUMBER.fullmatch(v) if isinstance(v, str) else None
    if m is None or (m.group(1) is not None and abs(int(m.group(1))) > 1000) or re.fullmatch(r".*/0+", v):
        raise ValueError(f"not an exact number: {v!r}")
    return Fraction(v)


def unique_keys(pairs):
    out = {}
    for k, v in pairs:
        if k in out:
            raise ValueError(f"duplicate key {k!r}")
        out[k] = v
    return out


BRACKETS = re.compile(r'\\.|["\[\]{}]', re.S)


def depth(v):
    """Nesting as brehm_cpp's reader counts it: every value is one level, the top one is level 1."""
    if isinstance(v, (dict, list)):
        return 1 + max(map(depth, v.values() if isinstance(v, dict) else v), default=0)
    return 1


def load(path):
    """Read a file with the reading rules of SPEC §3, including the 64-level nesting limit. Brackets are
    counted before parsing: with the high recursion limit set above, Python 3.11's JSON decoder would
    otherwise overflow the C stack on a deeply nested file."""
    with open(path, encoding="utf-8") as fh:
        text = fh.read()
    level, in_string = 0, False
    for m in BRACKETS.finditer(text):
        t = m.group()
        if t == '"':
            in_string = not in_string
        elif not in_string and t in "[{":
            level += 1
            if level > 64:
                raise ValueError("nesting deeper than 64 levels")
        elif not in_string and t in "]}":
            level -= 1
    data = json.loads(text, parse_float=num, parse_int=num, parse_constant=num, object_pairs_hook=unique_keys)
    if depth(data) > 64:
        raise ValueError("nesting deeper than 64 levels")
    return data


def pt(v):
    return (num(v[0]), num(v[1]))


def dot(u, v):
    return u[0] * v[0] + u[1] * v[1]


def sub(u, v):
    return (u[0] - v[0], u[1] - v[1])


def mat_vec(M, v):
    return (M[0][0] * v[0] + M[0][1] * v[1], M[1][0] * v[0] + M[1][1] * v[1])


def isometry(G):
    """Read the linear part of the starting isometry, checking G^T G = I exactly."""
    if not isinstance(G, (list, tuple)) or len(G) != 2 or any(
            not isinstance(row, (list, tuple)) or len(row) != 2 for row in G):
        raise ValueError("g must be a 2-by-2 matrix")
    G = tuple(tuple(num(x) for x in row) for row in G)
    if any(sum(G[k][i] * G[k][j] for k in range(2)) != int(i == j)
           for i in range(2) for j in range(2)):
        raise ValueError("g must be exactly orthogonal")
    return G


class Brehm:
    """f_1, ..., f_n of Brehm's construction for the pairs in insertion order."""

    def __init__(self, sources, targets, G, L, max_steps=400):
        self.a, self.b, self.G, self.L = sources, targets, isometry(G), L
        if L is not None and (L <= 0 or L * L < 2 * max(dot(sub(p, sources[0]), sub(p, sources[0])) for p in sources)):
            # Brehm's (II) and U in C (pp. 189-191) need d(a_i, a_1) <= r = L / sqrt(2)
            raise ValueError("L is not admissible: Brehm's construction needs L^2 >= 2 |a_i - a_1|^2")
        self.max_steps = max_steps
        self.memo = {}
        self.queries = 0           # evaluations of f_1
        self.max_offset = Fraction(0)  # largest |x - a_1|_inf at which f_1 was evaluated
        self.folded = False        # whether some f_1 query lay outside the inner square

    def f1(self, x):
        self.queries += 1
        a, L = self.a[0], self.L
        d = sub(x, a)
        self.max_offset = max(self.max_offset, abs(d[0]), abs(d[1]))
        s, phi = [], []
        for t in d:
            if L is not None and t > L:
                s.append(-1), phi.append(2 * L - t)
            elif L is not None and t < -L:
                s.append(-1), phi.append(-2 * L - t)
            else:
                s.append(1), phi.append(t)
        self.folded |= -1 in s
        G = self.G
        v = mat_vec(G, phi)
        value = (self.b[0][0] + v[0], self.b[0][1] + v[1])
        M = ((G[0][0] * s[0], G[0][1] * s[1]), (G[1][0] * s[0], G[1][1] * s[1]))
        Mx = mat_vec(M, x)
        return value, (M, (value[0] - Mx[0], value[1] - Mx[1]))

    def f(self, k, x):
        """(f_k(x), a germ near x); k is 1-based."""
        key = (k, x)
        if key in self.memo:
            return self.memo[key]
        if k == 1:
            out = self.f1(x)
        else:
            v, germ = self.f(k - 1, x)
            a, y = self.a[k - 1], self.b[k - 1]
            d = sub(x, a)
            if dot(sub(v, y), sub(v, y)) <= dot(d, d):
                out = (v, germ)  # x is not in U
            elif d == (0, 0):  # x = a: f_{i+1}(a) = b; borrow a germ (a search hint only) from a nearby point
                try:
                    hint = self.f(k, (a[0] + Fraction(1, 2**30), a[1] + Fraction(1, 2**31)))[1]
                except RuntimeError:
                    hint = None
                out = (y, hint)
            else:
                t, fp, (M, c) = self.exit(k - 1, a, y, d, germ)
                value = (y[0] + (fp[0] - y[0]) / t, y[1] + (fp[1] - y[1]) / t)
                # germ of the cone: (M, c) composed with the reflection fixing the bisector of a, q
                q = mat_vec(((M[0][0], M[1][0]), (M[0][1], M[1][1])), sub(y, c))
                u = sub(q, a)
                uu = dot(u, u)
                cc = (dot(q, q) - dot(a, a)) / 2
                S = ((1 - 2 * u[0] * u[0] / uu, -2 * u[0] * u[1] / uu), (-2 * u[1] * u[0] / uu, 1 - 2 * u[1] * u[1] / uu))
                MS = tuple(tuple(sum(M[i][m] * S[m][j] for m in range(2)) for j in range(2)) for i in range(2))
                shift = mat_vec(M, (2 * cc * u[0] / uu, 2 * cc * u[1] / uu))
                out = (value, (MS, (c[0] + shift[0], c[1] + shift[1])))
        self.memo[key] = out
        return out

    def exit(self, k, a, y, d, germ):
        """t > 1 with |f_k(a + t d) - y| = t |d| (exact), f_k there, and the germ used last.

        With finite L the search starts from the bracket [1, t_C], t_C where the ray leaves Brehm's
        protected square C = a_1 + [-3L, 3L]^2: U is inside the open square (Brehm p. 191, 4)), so
        the discrepancy is <= 0 there. The step budget grows with the bit length of the numbers, which
        bounds how close the exit can lie to a break of f_k along the ray."""
        t_lo, g_lo, t_hi, g_hi = Fraction(1), germ, None, None
        dd = dot(d, d)
        misses = 0
        budget = self.max_steps + 8 * max(max(v.numerator.bit_length(), v.denominator.bit_length())
                                          for v in (*a, *d, *y))
        if self.L is not None:
            a1, L = self.a[0], self.L
            t_hi = min(((a1[j] + (3 * L if d[j] > 0 else -3 * L)) - a[j]) / d[j] for j in range(2) if d[j] != 0)
            p = (a[0] + t_hi * d[0], a[1] + t_hi * d[1])
            fp, g_hi = self.f(k, p)
            s = dot(sub(fp, y), sub(fp, y)) - t_hi * t_hi * dd
            if s == 0:
                return t_hi, fp, g_hi
            if s > 0:
                raise RuntimeError("U reaches the boundary of the protected square, against Brehm's 4)")
            budget += 8 * max(t_hi.numerator.bit_length(), t_hi.denominator.bit_length())
        for _ in range(budget):
            if misses and misses % 16 == 0 and t_hi is not None:  # the germs may be stale: refit one from exact values (costly)
                g_lo = self.fit(k, a, d, t_lo, (t_hi - t_lo) / 1024) or g_lo
            cand, used = None, g_lo
            for g in (g_lo, g_hi):
                z = self.zero(a, y, d, g)
                if z is not None and z > t_lo and (t_hi is None or z < t_hi):
                    cand, used = z, g
                    break
            if cand is None:
                misses += 1
                cand = 2 * t_lo if t_hi is None else (t_lo + t_hi) / 2
            p = (a[0] + cand * d[0], a[1] + cand * d[1])
            fp, gp = self.f(k, p)
            s = dot(sub(fp, y), sub(fp, y)) - cand * cand * dd
            if s == 0:
                return cand, fp, used
            if s > 0:
                t_lo, g_lo = cand, gp
            else:
                t_hi, g_hi = cand, gp
        raise RuntimeError("no exit point found within the step budget"
                           + ("; with --L inf this can mean the value depends on the folds of f_1" if self.L is None else ""))

    @staticmethod
    def zero(a, y, d, g):
        """The t at which the germ g's discrepancy |g(a + t d) - y|^2 - t^2 |d|^2 vanishes, if any."""
        if g is None:
            return None
        M, c = g
        q = mat_vec(((M[0][0], M[1][0]), (M[0][1], M[1][1])), sub(y, c))
        u = sub(q, a)
        ud = dot(u, d)
        return dot(u, u) / (2 * ud) if ud > 0 else None

    def fit(self, k, a, d, t, eps):
        """The affine map through f_k at a + t d and two points eps further along d and across it,
        if it is an isometry (then, generically, the germ of f_k just beyond a + t d), else None."""
        x0 = (a[0] + t * d[0], a[1] + t * d[1])
        e1, e2 = (eps * d[0], eps * d[1]), (-eps * d[1], eps * d[0])
        v0 = self.f(k, x0)[0]
        v1 = self.f(k, (x0[0] + e1[0], x0[1] + e1[1]))[0]
        v2 = self.f(k, (x0[0] + e2[0], x0[1] + e2[1]))[0]
        det = e1[0] * e2[1] - e1[1] * e2[0]
        w1, w2 = sub(v1, v0), sub(v2, v0)
        # M [e1 e2] = [w1 w2]  =>  M = [w1 w2] [e1 e2]^-1
        inv = ((e2[1] / det, -e2[0] / det), (-e1[1] / det, e1[0] / det))
        M = tuple(tuple(w1[i] * inv[0][j] + w2[i] * inv[1][j] for j in range(2)) for i in range(2))
        if M[0][0] ** 2 + M[1][0] ** 2 != 1 or M[0][1] ** 2 + M[1][1] ** 2 != 1 or M[0][0] * M[0][1] + M[1][0] * M[1][1] != 0:
            return None
        Mx = mat_vec(M, x0)
        return M, (v0[0] - Mx[0], v0[1] - Mx[1])


def default_L(domain, a1):
    return 2 * max([Fraction(1)] + [abs(p[0] - a1[0]) + abs(p[1] - a1[1]) for p in domain])


def choices(data, how):
    """(sources, targets, order, g, L) with the free choices resolved as brehm_cpp resolves them;
    `how` holds the optional "order", "g" and "L"."""
    domain = [pt(v) for v in data["domain"]]
    sources = [pt(v) for v in data["sources"]]
    targets = [pt(v) for v in data["targets"]]
    order = [num(i) for i in how.get("order", [str(i) for i in range(len(sources))])]
    if any(k.denominator != 1 for k in order) or sorted(order) != list(range(len(sources))):
        raise ValueError("order is not a permutation of the sources")
    order = [int(k) for k in order]
    G = isometry(how.get("g", [["1", "0"], ["0", "1"]]))
    L = num(how["L"]) if "L" in how else default_L(domain, sources[order[0]])
    return sources, targets, order, G, L, domain


def from_certificate(cert, L_override=None):
    sources, targets, order, G, L, domain = choices(cert, cert.get("construction", {}))
    a = [sources[i] for i in order]
    b = [targets[i] for i in order]
    return Brehm(a, b, G, None if L_override == "inf" else L), domain


def apply(piece, x):
    A = [pt(r) for r in piece["A"]]
    return (A[0][0] * x[0] + A[0][1] * x[1] + num(piece["b"][0]), A[1][0] * x[0] + A[1][1] * x[1] + num(piece["b"][1]))


def faithful(cert, L_override=None, points=0, seed=0, instance=None):
    """Compare the certificate with f_n at every cell vertex (and at `points` random points). With an
    instance, the pairs and the resolved choices must also be the instance's."""
    if instance is not None and choices(cert, cert.get("construction", {}))[:5] != choices(instance, instance)[:5]:
        return {"agree": False, "mismatches": [], "instance": "the certificate's pairs or choices are not the instance's"}
    br, domain = from_certificate(cert, L_override)
    n = len(br.a)
    pieces = [(p, [pt(v) for v in p["polygon"]]) for p in cert["pieces"]]
    mismatches, checked = [], 0
    for idx, (piece, poly) in enumerate(pieces):
        for v in poly:
            checked += 1
            if br.f(n, v)[0] != apply(piece, v):
                mismatches.append({"cell": idx, "vertex": [str(c) for c in v]})
    rng = random.Random(seed)
    xs = [p[0] for p in domain]
    ys = [p[1] for p in domain]
    sampled = 0
    while sampled < points:
        x = (Fraction(rng.randint(0, 10**6), 10**6) * (max(xs) - min(xs)) + min(xs),
             Fraction(rng.randint(0, 10**6), 10**6) * (max(ys) - min(ys)) + min(ys))
        for idx, (piece, poly) in enumerate(pieces):
            if all((poly[i - 1][0] - x[0]) * (poly[i][1] - x[1]) - (poly[i - 1][1] - x[1]) * (poly[i][0] - x[0]) >= 0
                   for i in range(len(poly))) or all(
                    (poly[i - 1][0] - x[0]) * (poly[i][1] - x[1]) - (poly[i - 1][1] - x[1]) * (poly[i][0] - x[0]) <= 0
                    for i in range(len(poly))):
                sampled += 1
                if br.f(n, x)[0] != apply(piece, x):
                    mismatches.append({"cell": idx, "point": [str(c) for c in x]})
                break
    images = all(br.f(n, s)[0] == t for s, t in zip(br.a, br.b))
    return {"vertices_checked": checked, "random_points": sampled, "mismatches": mismatches[:25],
            "sources_hit_targets": images, "L": "inf" if br.L is None else str(br.L),
            "f1_queries": br.queries, "max_offset": str(br.max_offset), "fold_branch_used": br.folded,
            "agree": not mismatches and images}


def main():
    args = sys.argv[1:]
    if not args:
        print(__doc__, file=sys.stderr)
        return 2
    opts = dict(zip(args[1::2], args[2::2]))
    if len(args) % 2 == 0 or set(opts) - {"--instance", "--L", "--points", "--seed"}:
        print(__doc__, file=sys.stderr)
        return 2
    try:
        cert = load(args[0])
        instance = load(opts["--instance"]) if "--instance" in opts else None
    except (OSError, ValueError, TypeError, RecursionError) as e:
        print(f"invalid input: {e!r}", file=sys.stderr)
        return 3
    try:
        rep = faithful(cert, opts.get("--L"), int(opts.get("--points", 0)), int(opts.get("--seed", 0)), instance)
    except (ValueError, KeyError, TypeError, AttributeError, IndexError, ZeroDivisionError) as e:
        print(f"invalid input: {e!r}", file=sys.stderr)
        return 3
    except MemoryError:
        print("out of memory", file=sys.stderr)
        return 5
    except Exception as e:  # noqa: BLE001 - e.g. no exit point within the step budget: undecided, not a verdict
        print(f"undecided: {e!r}", file=sys.stderr)
        return 4
    print(json.dumps(rep))
    return 0 if rep["agree"] else 1


if __name__ == "__main__":
    sys.set_int_max_str_digits(0)
    sys.exit(main())
