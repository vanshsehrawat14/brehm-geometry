"""Reference checker for brehm-rational-v1 certificates (docs/SPEC.md, conditions R1-R6).

Standard library only; shares no code with the constructors or with the other two checkers
(python/core/exact_certify.py: pairwise overlap + area sum; check.cpp: boundary cancellation).
Coverage is checked by a third method, a vertical slab decomposition: between consecutive vertex
x-coordinates every cell is a trapezoid bounded by two lines, and the cells must stack from the
domain's lower line to its upper line with equal lines at every contact. Continuity is checked on
every pair of intersecting cells at all extreme points of the intersection (not only on shared
edges).

    python ref_check.py CERT.json [--instance INSTANCE.json]

With --instance it also checks that the certificate's domain and prescribed pairs are those of the
instance. Prints a JSON report; exit code 0 if every condition holds, 1 if not, 2 on bad usage,
3 on unreadable input (including duplicate JSON keys and numbers outside the grammar of SPEC §3),
4 on an internal error, 5 when out of memory. All arithmetic is exact (fractions.Fraction); JSON
numbers and strings are read exactly.
"""
import json
import re
import sys
from fractions import Fraction

NUMBER = re.compile(r"[+-]?(?:\d+/\d+|(?:\d+\.?\d*|\.\d+)(?:[eE]([+-]?\d+))?)")


def num(v):
    """An exact number: a Fraction, or a string that is an integer, a decimal with |exponent| <= 1000,
    or p/q with q > 0. The same grammar as brehm_cpp (json.hpp)."""
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


def point(v):
    if not isinstance(v, list) or len(v) != 2:
        raise ValueError(f"not a point: {v!r}")
    return (num(v[0]), num(v[1]))


BRACKETS = re.compile(r'\\.|["\[\]{}]', re.S)


def depth(v):
    """Nesting as brehm_cpp's reader counts it: every value is one level, the top one is level 1."""
    if isinstance(v, (dict, list)):
        return 1 + max(map(depth, v.values() if isinstance(v, dict) else v), default=0)
    return 1


def load(path):
    """Read a file with the reading rules of SPEC §3, including the 64-level nesting limit. Brackets are
    counted before parsing, so a deeply nested file never reaches the recursive parser."""
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


def turn(o, a, b):
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def convex_ccw(poly):
    """The polygon as a counter-clockwise list with collinear vertices removed, or None if it is
    not a convex polygon with positive area listed in boundary order."""
    if len(poly) >= 2 and poly[0] == poly[-1]:
        poly = poly[:-1]
    if len(poly) < 3 or len(set(poly)) != len(poly):
        return None
    area = sum(turn((0, 0), poly[i], poly[(i + 1) % len(poly)]) for i in range(len(poly)))
    if area == 0:
        return None
    if area < 0:
        poly = poly[::-1]
    n = len(poly)
    for i in range(n):  # every vertex weakly left of every edge, and the boundary winds once
        for w in poly:
            if turn(poly[i], poly[(i + 1) % n], w) < 0:
                return None
    out = [poly[i] for i in range(n) if turn(poly[i - 1], poly[i], poly[(i + 1) % n]) != 0]
    return out if len(out) >= 3 else None


def inside(ccw, p):
    return all(turn(ccw[i], ccw[(i + 1) % len(ccw)], p) >= 0 for i in range(len(ccw)))


def apply(cell, p):
    (a, b), (c, d) = cell["A"]
    return (a * p[0] + b * p[1] + cell["b"][0], c * p[0] + d * p[1] + cell["b"][1])


def slab_lines(ccw, xl, xr):
    """The (lower, upper) lines y = m x + k of the convex polygon over the open slab (xl, xr)."""
    lines = []
    for i in range(len(ccw)):
        p, q = ccw[i], ccw[(i + 1) % len(ccw)]
        if p[0] != q[0] and min(p[0], q[0]) <= xl and max(p[0], q[0]) >= xr:
            m = (q[1] - p[1]) / (q[0] - p[0])
            lines.append((m, p[1] - m * p[0]))
    assert len(lines) == 2, "a convex polygon spanning a slab has one lower and one upper edge there"
    xm = (xl + xr) / 2
    lines.sort(key=lambda l: l[0] * xm + l[1])
    return lines[0], lines[1]


def check_tiling(domain, polys, fail):
    """R3 by slabs: inside every open vertical slab between consecutive vertex x-coordinates, the
    cells meeting it must stack exactly from the domain's lower line to its upper line."""
    xs = sorted({p[0] for poly in polys for p in poly} | {p[0] for p in domain})
    span = [(min(p[0] for p in poly), max(p[0] for p in poly)) for poly in polys]
    order = sorted(range(len(polys)), key=lambda i: span[i][0])
    dlo, dhi = min(p[0] for p in domain), max(p[0] for p in domain)
    active, nxt, slabs = [], 0, 0
    for xl, xr in zip(xs, xs[1:]):
        while nxt < len(order) and span[order[nxt]][0] <= xl:
            active.append(order[nxt])
            nxt += 1
        active = [i for i in active if span[i][1] >= xr]
        slabs += 1
        if not dlo <= xl < dhi:
            if active:
                fail(f"cell {active[0]} lies outside the domain near x = {xl}")
            continue
        xm = (xl + xr) / 2
        stack = sorted(((slab_lines(polys[i], xl, xr), i) for i in active),
                       key=lambda t: t[0][0][0] * xm + t[0][0][1])
        low, high = slab_lines(domain, xl, xr)
        cur = low
        for (lo, hi), i in stack:
            if lo != cur:
                fail(f"gap or overlap below cell {i} in the slab ({xl}, {xr})")
                break
            cur = hi
        else:
            if cur != high:
                fail(f"the cells do not reach the top of the domain in the slab ({xl}, {xr})")
    return slabs


def segment_points(p1, p2, q1, q2):
    d = turn((0, 0), (p2[0] - p1[0], p2[1] - p1[1]), (q2[0] - q1[0], q2[1] - q1[1]))
    if d == 0:
        return []  # parallel: the extreme points of a collinear overlap are endpoints
    w = (q1[0] - p1[0], q1[1] - p1[1])
    t = turn((0, 0), w, (q2[0] - q1[0], q2[1] - q1[1])) / d
    s = turn((0, 0), w, (p2[0] - p1[0], p2[1] - p1[1])) / d
    if 0 <= t <= 1 and 0 <= s <= 1:
        return [(p1[0] + t * (p2[0] - p1[0]), p1[1] + t * (p2[1] - p1[1]))]
    return []


def check_continuity(polys, cells, fail):
    """R5: for every pair of cells, the two maps agree on their intersection. The intersection of
    two convex polygons is convex and its extreme points are vertices of one polygon lying in the
    other or crossings of two edges; two affine maps that agree there agree on the whole set."""
    box = [(min(p[0] for p in q), max(p[0] for p in q), min(p[1] for p in q), max(p[1] for p in q)) for q in polys]
    order = sorted(range(len(polys)), key=lambda i: box[i][0])
    pairs = 0
    for a, i in enumerate(order):
        for j in order[a + 1:]:
            if box[j][0] > box[i][1]:
                break
            if box[j][3] < box[i][2] or box[j][2] > box[i][3]:
                continue
            P, Q = polys[i], polys[j]
            pts = [v for v in P if inside(Q, v)] + [v for v in Q if inside(P, v)]
            for k in range(len(P)):
                for m in range(len(Q)):
                    pts += segment_points(P[k], P[(k + 1) % len(P)], Q[m], Q[(m + 1) % len(Q)])
            if pts:
                pairs += 1
            for v in pts:
                if apply(cells[i], v) != apply(cells[j], v):
                    fail(f"cells {i} and {j} disagree at {tuple(map(str, v))}")
                    break
    return pairs


def same_polygon(p, q):
    return p is not None and q is not None and len(p) == len(q) and any(
        p == q[k:] + q[:k] for k in range(len(q)))


def check(cert, instance=None):
    failures = []

    def fail(msg):
        failures.append(msg)

    report = {"format": cert.get("format")}
    if cert.get("format") != "brehm-rational-v1":
        fail("unsupported certificate format")
        return report, failures
    domain = convex_ccw([point(v) for v in cert["domain"]])  # R1
    if domain is None:
        fail("R1: the domain is not a convex polygon with positive area")
        return report, failures
    sources = [point(v) for v in cert["sources"]]
    targets = [point(v) for v in cert["targets"]]
    if len(sources) != len(targets):
        fail("source/target length mismatch")
    cells, polys = [], []
    for k, raw in enumerate(cert["pieces"]):
        A = [point(r) for r in raw["A"]]
        if len(A) != 2:
            fail(f"R4: cell {k}: the linear part is not 2-by-2")
            continue
        cell = {"A": A, "b": point(raw["b"])}
        (a, b), (c, d) = A
        if a * a + c * c != 1 or b * b + d * d != 1 or a * b + c * d != 0:
            fail(f"R4: cell {k}: the linear part is not orthogonal")
        poly = convex_ccw([point(v) for v in raw["polygon"]])  # R2
        if poly is None:
            fail(f"R2: cell {k} is not a convex polygon with positive area")
            continue
        cells.append(cell)
        polys.append(poly)
    report["cells"] = len(polys)
    if failures:
        return report, failures
    report["slabs"] = check_tiling(domain, polys, fail)  # R3
    if failures:
        return report, failures
    report["intersecting_pairs"] = check_continuity(polys, cells, fail)  # R5
    for k, (s, t) in enumerate(zip(sources, targets)):  # R6
        if not inside(domain, s):
            fail(f"R6: source {k} is outside the domain")
        hits = [i for i, p in enumerate(polys) if inside(p, s)]
        if not hits:
            fail(f"R6: source {k} is not covered")
        for i in hits:
            if apply(cells[i], s) != t:
                fail(f"R6: cell {i} sends source {k} to {tuple(map(str, apply(cells[i], s)))}, not to its target")
    if instance is not None:
        want = convex_ccw([point(v) for v in instance["domain"]])
        if not same_polygon(domain, want):
            fail("the certificate's domain is not the instance's domain")
        pairs = sorted(zip([point(v) for v in instance["sources"]], [point(v) for v in instance["targets"]]))
        if pairs != sorted(zip(sources, targets)):
            fail("the certificate's prescribed pairs are not the instance's")
    return report, failures


def main():
    args = sys.argv[1:]
    if not args or len(args) not in (1, 3) or (len(args) == 3 and args[1] != "--instance"):
        print(__doc__, file=sys.stderr)
        return 2
    try:
        cert = load(args[0])
        instance = load(args[2]) if len(args) == 3 else None
        report, failures = check(cert, instance)
    except (OSError, ValueError, KeyError, TypeError, AttributeError, IndexError, RecursionError) as e:
        print(f"invalid input: {e!r}", file=sys.stderr)
        return 3
    except MemoryError:
        print("out of memory", file=sys.stderr)
        return 5
    except Exception as e:  # noqa: BLE001 - any other failure is ours, not a verdict on the input
        print(f"internal error: {e!r}", file=sys.stderr)
        return 4
    report["all_pass"] = not failures
    report["failures"] = failures[:25]
    print(json.dumps(report))
    return 0 if not failures else 1


if __name__ == "__main__":
    sys.set_int_max_str_digits(0)  # long exact numbers are legitimate; brehm_cpp has no limit either
    sys.exit(main())
