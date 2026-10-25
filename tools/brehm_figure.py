"""Draw a brehm-rational-v1 certificate: the paper cut into its cells, and the images of the cells.

    python brehm_figure.py CERT.json OUT.svg [--exe BREHM_CPP | --unverified]
                           [--labels | --no-labels] [--highlight I,J,...] [--zoom I,J,...] [--no-points]

The certificate is first checked by `brehm_cpp check` (the executable given by --exe, else the one on
PATH); a rejected certificate gives no figure. --unverified skips that check, and the figure says so.

Left panel: the paper D and the certificate's cells, tinted by the orientation of each cell's
isometry (det A = +1 or -1), and the sources a_k. Where two cells meet along an edge, the edge is a
crease if their maps differ (they then differ by the reflection in that edge) and an auxiliary edge if
the maps are equal; creases, auxiliary edges and the paper's boundary are drawn differently. Right
panel, at the same scale: the image A x + b of every cell, the images of the same edges, and the
targets b_k. Every image cell gets the same translucent fill, so a point looks darker the more cells
cover it, whatever the drawing order: no layer order is drawn, because a certificate contains none.
Mountain and valley are not assigned. In a browser, hovering over a cell shows its index.

All geometry is exact (fractions.Fraction; numbers are read with the rules of docs/SPEC.md §3). Each
coordinate is rounded once, to 0.001 px, when it is written. The figure illustrates the certificate;
the certificate is the mathematical object. With brehm_cpp's check, the edge structure computed here
must agree with the checker's own counts, or no figure is written. Standard library only.

Exit codes: 0 figure written; 1 certificate rejected (by brehm_cpp check, or cells whose edges do not
fit together); 2 usage, or brehm_cpp not found; 3 unreadable input; 4 internal error; 5 out of memory.
"""
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
import time
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path
from xml.sax.saxutils import escape

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ref_check import load, num, point  # noqa: E402  (the reading rules of docs/SPEC.md §3)

PANEL, GAP, LEFT, TOP = 460, 40, 20, 70        # px: panel side, gap between panels, margins
WIDTH = LEFT + 2 * PANEL + GAP + LEFT
FRONT, BACK = "#f4ecd2", "#c9dbee"              # paper tints: det A = +1, det A = -1
INK, ALPHA = (27, 42, 65), Fraction(1, 5)       # the one fill of every image cell, and its opacity
STYLE = {"aux": 'stroke="#808080" stroke-width="0.6" stroke-dasharray="3 2"',
         "edge": 'stroke="#000000" stroke-width="1.3"',
         "crease": 'stroke="#c0392b" stroke-width="1.1"'}
HIGHLIGHT = ("#e66101", "#5e3c99", "#1a9641", "#d01c8b", "#0571b0", "#b8860b")
AUTO_LABELS = 60                                 # cell labels by default up to this many cells


class Rejected(Exception):
    """The cells cannot be drawn as a map: their edges do not fit together."""


def read(raw):
    if raw.get("format") != "brehm-rational-v1":
        raise ValueError("unsupported certificate format")
    domain = [point(v) for v in raw["domain"]]
    sources = [point(v) for v in raw["sources"]]
    targets = [point(v) for v in raw["targets"]]
    if len(sources) != len(targets):
        raise ValueError("sources and targets differ in length")
    cells = []
    for k, c in enumerate(raw["pieces"]):
        rows = [point(r) for r in c["A"]]
        if len(rows) != 2:
            raise ValueError(f"cell {k}: the linear part is not 2-by-2")
        cells.append(([point(v) for v in c["polygon"]], (rows[0], rows[1]), point(c["b"])))
    return domain, sources, targets, cells


def apply(cell, p):
    _, ((a, b), (c, d)), t = cell
    return (a * p[0] + b * p[1] + t[0], c * p[0] + d * p[1] + t[1])


def orientation(poly):
    s = sum(p[0] * q[1] - q[0] * p[1] for p, q in zip(poly, poly[1:] + poly[:1]))
    return (s > 0) - (s < 0)


def edge_records(poly, owner, sense):
    """Each edge on its exact supporting line, as brehm_cpp's check groups them: kind 1 is the line
    x + k1 y = k2 with parameter y, kind 0 the horizontal line y = k2 with parameter x. The sign is +1
    if the edge runs the way it would on the counter-clockwise boundary (sense = orientation)."""
    for p, q in zip(poly, poly[1:] + poly[:1]):
        if p == q:
            raise Rejected(f"{'the paper' if owner < 0 else f'cell {owner}'} repeats a vertex")
        dx, dy = q[0] - p[0], q[1] - p[1]
        if dy:
            k1 = -dx / dy
            yield (1, k1, p[0] + k1 * p[1]), min(p[1], q[1]), max(p[1], q[1]), (1 if dy > 0 else -1) * sense, owner
        else:
            yield (0, Fraction(0), p[1]), min(p[0], q[0]), max(p[0], q[0]), (1 if dx > 0 else -1) * sense, owner


def at(key, s):
    return (key[2] - key[1] * s, s) if key[0] == 1 else (s, key[2])


def structure(domain, cells):
    """The elementary edge intervals: on every line, edges are split at all their endpoints. In a
    valid certificate each interval has exactly two edges of opposite sense: two cells (a crease or
    an auxiliary edge) or one cell and the paper's boundary. Returns [(kind, u, v, p0, p1)]."""
    groups = {}
    for i, cell in enumerate(cells):
        o = orientation(cell[0])
        if len(cell[0]) < 3 or o == 0:
            raise Rejected(f"cell {i} is not a polygon with positive area")
        for key, *rec in edge_records(cell[0], i, o):
            groups.setdefault(key, []).append(rec)
    o = orientation(domain)
    if len(domain) < 3 or o == 0:
        raise Rejected("the paper is not a polygon with positive area")
    for key, *rec in edge_records(domain, -1, -o):
        groups.setdefault(key, []).append(rec)
    out = []
    for key in sorted(groups):
        recs = groups[key]
        cuts = sorted({c for lo, hi, _, _ in recs for c in (lo, hi)})
        by_lo = sorted(range(len(recs)), key=lambda r: recs[r][0])
        by_hi = sorted(range(len(recs)), key=lambda r: recs[r][1])
        active, ia, ib = set(), 0, 0
        for c0, c1 in zip(cuts, cuts[1:]):
            while ia < len(recs) and recs[by_lo[ia]][0] <= c0:
                active.add(by_lo[ia])
                ia += 1
            while ib < len(recs) and recs[by_hi[ib]][1] <= c0:
                active.discard(by_hi[ib])
                ib += 1
            if not active:
                continue
            p0, p1 = at(key, c0), at(key, c1)
            owners = sorted(recs[r][3] for r in active)
            if len(owners) != 2 or sum(recs[r][2] for r in active) != 0 or owners[0] == owners[1]:
                raise Rejected(f"the cells do not fit together near ({p0[0]}, {p0[1]}): "
                               f"{len(owners)} edges there")
            u, v = owners
            if u < 0:
                out.append(("edge", v, None, p0, p1))
                continue
            if apply(cells[u], p0) != apply(cells[v], p0) or apply(cells[u], p1) != apply(cells[v], p1):
                raise Rejected(f"cells {u} and {v} map their common edge differently")
            same = cells[u][1:] == cells[v][1:]
            out.append(("aux" if same else "crease", u, v, p0, p1))
    return out


def faces(n, intervals):
    """Components of cells glued across auxiliary edges (brehm_cpp check's isometric regions)."""
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for kind, u, v, _, _ in intervals:
        if kind == "aux":
            parent[find(u)] = find(v)
    return sum(find(i) == i for i in range(n))


def det(cell):
    (a, b), (c, d) = cell[1]
    return a * d - b * c


# --- drawing -----------------------------------------------------------------------------------

def fmt(q):
    """q rounded to 0.001 (ties to even), as a short decimal; exact, no floating point."""
    n = round(q * 1000)
    i, f = divmod(abs(n), 1000)
    return ("-" if n < 0 else "") + str(i) + (f".{f:03d}".rstrip("0") if f else "")


def approx(q, digits=4):
    with localcontext() as ctx:
        ctx.prec = digits
        d = Decimal(q.numerator) / Decimal(q.denominator)
    return format(d.normalize(), "f" if -5 <= d.adjusted() <= 6 else "e").replace("-", "−")


def show(q, width=10):
    s = str(q).replace("-", "−")
    return s if len(s) <= width else "≈" + approx(q)


def bbox(points):
    xs, ys = [p[0] for p in points], [p[1] for p in points]
    return min(xs), min(ys), max(xs), max(ys)


class View:
    """Exact map from the plane to one panel: X = ox + s x, Y = oy - s y, the box centred."""

    def __init__(self, box, s, x0, y0):
        xmin, ymin, xmax, ymax = box
        self.s = s
        self.ox = x0 + (PANEL - (xmax - xmin) * s) / 2 - xmin * s
        self.oy = y0 + (PANEL - (ymax - ymin) * s) / 2 + ymax * s

    def __call__(self, p):
        return fmt(self.ox + self.s * p[0]), fmt(self.oy - self.s * p[1])

    def pts(self, poly):
        return " ".join("%s,%s" % self(p) for p in poly)


def nice_length(limit):
    """The largest of 1, 2, 5 times a power of ten that is at most limit (> 0)."""
    k = len(str(limit.numerator)) - len(str(limit.denominator))
    while Fraction(10) ** k > limit:
        k -= 1
    while Fraction(10) ** (k + 1) <= limit:
        k += 1
    m = max(m for m in (1, 2, 5) if m * Fraction(10) ** k <= limit)
    return m * Fraction(10) ** k


def shade(k):
    """The colour of k stacked image fills over white, as a viewer composites them."""
    c = [round(v + (1 - ALPHA) ** k * (255 - v)) for v in INK]
    return "#%02x%02x%02x" % tuple(c)


def ranges(idx):
    idx, out, i = sorted(idx), [], 0
    while i < len(idx):
        j = i
        while j + 1 < len(idx) and idx[j + 1] == idx[j] + 1:
            j += 1
        out.append(str(idx[i]) if j == i else f"{idx[i]}–{idx[j]}" if j > i + 1 else f"{idx[i]},{idx[j]}")
        i = j + 1
    return ",".join(out)


def grouped(points):
    """Equal points drawn once, labelled with all their indices."""
    at_point = {}
    for k, p in enumerate(points):
        at_point.setdefault(p, []).append(k)
    return sorted(at_point.items(), key=lambda kv: kv[1][0])


def cells_named(ks):
    return ("cell " if len(ks) == 1 else "cells ") + ", ".join(map(str, ks))


def describe_choices(how):
    """The recorded construction choices, for the caption only (the checkers ignore them)."""
    try:
        order = [str(num(k)) for k in how["order"]]
        g = [point(row) for row in how["g"]]
        L = num(how["L"])
    except (KeyError, TypeError, ValueError):
        return "no construction record"
    gs = "id" if g == [(1, 0), (0, 1)] else "[" + "; ".join(" ".join(show(x, 8) for x in row) for row in g) + "]"
    order_s = ",".join(order) if len(order) <= 16 else ",".join(order[:12]) + ",…"
    return f"order {order_s}; g = {gs}; L = {show(L, 12)}"


def draw(name, digest, domain, sources, targets, cells, intervals, how, verified, labels, highlight, points, zoom):
    m = len(cells)
    images = [[apply(c, p) for p in c[0]] for c in cells]
    if zoom:  # fit both panels to these cells; everything else is still drawn, clipped to the panels
        paper_box = bbox([p for k in zoom for p in cells[k][0]])
        image_box = bbox([p for k in zoom for p in images[k]])
    else:
        paper_box = bbox(domain + (sources if points else []))
        image_box = bbox([p for poly in images for p in poly] + (targets if points else []))
    span = max(paper_box[2] - paper_box[0], paper_box[3] - paper_box[1],
               image_box[2] - image_box[0], image_box[3] - image_box[1])
    s = Fraction(PANEL) / span
    x_left, x_right = LEFT, LEFT + PANEL + GAP
    left, right = View(paper_box, s, x_left, TOP), View(image_box, s, x_right, TOP)
    count = {k: sum(1 for t in intervals if t[0] == k) for k in ("crease", "aux", "edge")}
    nfaces = faces(m, intervals)
    height = TOP + PANEL + 184
    ink = "#%02x%02x%02x" % INK
    svg = ['<?xml version="1.0" encoding="UTF-8"?>',
           f'<svg xmlns="http://www.w3.org/2000/svg" version="1.1" width="{WIDTH}" height="{height}" '
           f'viewBox="0 0 {WIDTH} {height}" font-family="sans-serif" font-size="12">',
           f"<title>{escape(name)}</title>",
           "<desc>" + escape(json.dumps({
               "tool": "brehm_figure.py", "certificate": name, "sha256": digest, "cells": m,
               "verified": verified, "zoom": zoom, "scale_px_per_unit": str(s),
               "paper_transform": [str(left.ox), str(left.oy)], "image_transform": [str(right.ox), str(right.oy)],
               "rounding": "every coordinate rounded to 0.001 px"})) + "</desc>",
           # reversed image cells are hatched in the fill colour; layers of one colour composite in any
           # order to the same result, so the drawing order cannot show
           f'<defs><pattern id="reversed" patternUnits="userSpaceOnUse" width="5" height="5" '
           f'patternTransform="rotate(45)"><path d="M0 0v5" stroke="{ink}" stroke-width="1.4" '
           f'stroke-opacity="0.45"/></pattern>'
           f'<clipPath id="paper-panel"><rect x="{x_left}" y="{TOP}" width="{PANEL}" height="{PANEL}"/></clipPath>'
           f'<clipPath id="image-panel"><rect x="{x_right}" y="{TOP}" width="{PANEL}" height="{PANEL}"/></clipPath>'
           "</defs>",
           f'<rect width="{WIDTH}" height="{height}" fill="#ffffff"/>']
    title = f"{name}: {len(sources)} pairs, {m} cells, {nfaces} faces (regions with one isometry)"
    svg.append(f'<text x="{LEFT}" y="22" font-size="15" font-weight="bold">{escape(title)}</text>')
    status = ("certificate accepted by brehm_cpp check (R1–R6)" if verified else
              "NOT VERIFIED: drawn without running brehm_cpp check")
    zoomed = f" · zoomed to {cells_named(zoom)}" if zoom else ""
    warn = "" if verified else ' fill="#b00000"'
    svg.append(f'<text x="{LEFT}" y="42" font-size="12"{warn}>'
               f"{escape(describe_choices(how))} · {escape(status)}{escape(zoomed)}</text>")
    for x0, caption in ((x_left, "Paper D cut into the certificate's cells"),
                        (x_right, "Image F(D) at the same scale: all cells translucent, no layer order")):
        svg.append(f'<text x="{x0}" y="{TOP - 8}" font-size="12">{escape(caption)}</text>')
        svg.append(f'<rect x="{x0}" y="{TOP}" width="{PANEL}" height="{PANEL}" fill="none" stroke="#d0d0d0"/>')

    def cell_title(k):
        return f"<title>cell {k}{'' if det(cells[k]) > 0 else ' (reversed)'}</title>"

    def panel(view, image):
        out = []
        if not image:  # the paper: cells tinted by orientation
            for front, tint in ((True, FRONT), (False, BACK)):
                out.append(f'<g fill="{tint}" stroke="none">')
                out += [f'<polygon points="{view.pts(c[0])}">{cell_title(k)}</polygon>'
                        for k, c in enumerate(cells) if (det(c) > 0) == front]
                out.append("</g>")
        else:  # the image: one translucent fill for every cell, reversed cells also hatched
            out.append(f'<g fill="{ink}" fill-opacity="{float(ALPHA)}" stroke="none">')
            out += [f'<polygon points="{view.pts(poly)}">{cell_title(k)}</polygon>' for k, poly in enumerate(images)]
            out.append('</g><g fill="url(#reversed)" stroke="none">')
            out += [f'<polygon points="{view.pts(poly)}"/>' for k, poly in enumerate(images) if det(cells[k]) < 0]
            out.append("</g>")
        for kind in ("aux", "edge", "crease"):  # a fixed precedence where drawn lines coincide
            d = []
            for t in intervals:
                if t[0] == kind:
                    p0, p1 = (apply(cells[t[1]], t[3]), apply(cells[t[1]], t[4])) if image else (t[3], t[4])
                    d.append("M%s %sL%s %s" % (view(p0) + view(p1)))
            if d:
                out.append(f'<path fill="none" stroke-linecap="round" {STYLE[kind]} d="{"".join(d)}"/>')
        # outlines: a hairline on each boundary (it shows even a cell too thin to draw), then a band inside
        # each cell, clipped to it, so two outlined cells that share an edge each keep their colour on
        # their own side
        outlined = [(HIGHLIGHT[n % len(HIGHLIGHT)], view.pts(images[k] if image else cells[k][0]))
                    for n, k in enumerate(highlight)]
        out += [f'<polygon points="{p}" fill="none" stroke="{c}" stroke-width="1"/>' for c, p in outlined]
        for n, (c, p) in enumerate(outlined):
            clip = f"{'image' if image else 'paper'}-outline-{n}"
            out.append(f'<clipPath id="{clip}"><polygon points="{p}"/></clipPath>'
                       f'<polygon points="{p}" fill="none" stroke="{c}" stroke-width="5" clip-path="url(#{clip})"/>')
        if points:
            for p, ks in grouped(targets if image else sources):
                x, y = view(p)
                out.append(f'<circle cx="{x}" cy="{y}" r="2.6" fill="#000000"/>')
                out.append(f'<text x="{x}" y="{y}" dx="4" dy="-4" font-size="11">{ranges(ks)}</text>')
        if labels:
            at_label = {}
            for k, c in enumerate(cells):
                mid = (sum(p[0] for p in c[0]) / len(c[0]), sum(p[1] for p in c[0]) / len(c[0]))
                at_label.setdefault(view(apply(c, mid) if image else mid), []).append(k)
            for (x, y), ks in sorted(at_label.items(), key=lambda kv: kv[1][0]):
                out.append(f'<text x="{x}" y="{y}" dy="0.35em" text-anchor="middle" font-size="9" '
                           f'fill="#404040">{ranges(ks)}</text>')
        return out

    svg.append('<g id="paper" clip-path="url(#paper-panel)">')
    svg += panel(left, False)
    svg.append('</g>\n<g id="image" clip-path="url(#image-panel)">')
    svg += panel(right, True)
    svg.append("</g>")
    # scale bars and coordinate ranges
    bar = nice_length(Fraction(PANEL, 4) / s)
    for x0, box in ((x_left, paper_box), (x_right, image_box)):
        y = TOP + PANEL + 18
        svg.append(f'<path d="M{x0} {y}h{fmt(bar * s)}M{x0} {y - 4}v8M{fmt(x0 + bar * s)} {y - 4}v8" '
                   f'stroke="#000000" stroke-width="1.2" fill="none"/>')
        svg.append(f'<text x="{fmt(x0 + bar * s + 6)}" y="{y + 4}" font-size="11">{escape(show(bar))}</text>')
        rng = (f"x from {show(box[0])} to {show(box[2])}, y from {show(box[1])} to {show(box[3])}")
        svg.append(f'<text x="{x0 + PANEL}" y="{y + 4}" font-size="11" text-anchor="end">{escape(rng)}</text>')
    # legend
    y = TOP + PANEL + 46
    items = [("rect", FRONT, "isometry keeps orientation (det A = +1)"),
             ("rect", BACK, "isometry reverses orientation (det A = −1)"),
             ("line", "edge", "boundary of the paper and its image"),
             ("line", "crease", f"crease: the two maps differ ({count['crease']} pieces); M/V not assigned"),
             ("line", "aux", f"auxiliary edge: one map on both sides ({count['aux']} pieces)")]
    for n, (what, style, text) in enumerate(items):
        yy = y + 16 * n
        if what == "rect":
            svg.append(f'<rect x="{LEFT}" y="{yy - 9}" width="22" height="11" fill="{style}" stroke="#a0a0a0" stroke-width="0.5"/>')
        else:
            svg.append(f'<path d="M{LEFT} {yy - 4}h22" fill="none" {STYLE[style]}/>')
        svg.append(f'<text x="{LEFT + 30}" y="{yy}">{escape(text)}</text>')
    svg.append(f'<text x="{x_right}" y="{y}">Cells covering a point (shading only; it saturates):</text>')
    for n, k in enumerate((1, 2, 3, 4, 6)):
        svg.append(f'<rect x="{x_right + 46 * n}" y="{y + 7}" width="40" height="13" fill="{shade(k)}"/>')
        svg.append(f'<text x="{x_right + 46 * n + 20}" y="{y + 34}" text-anchor="middle" font-size="11">{k}</text>')
    svg.append(f'<rect x="{x_right}" y="{y + 43}" width="22" height="11" fill="url(#reversed)" stroke="#a0a0a0" stroke-width="0.5"/>')
    svg.append(f'<text x="{x_right + 30}" y="{y + 52}">hatched: image of a reversed cell (blue on the left)</text>')
    if highlight:
        marks = ", ".join(f'<tspan fill="{HIGHLIGHT[n % len(HIGHLIGHT)]}" font-weight="bold">{k}</tspan>'
                          for n, k in enumerate(highlight))
        svg.append(f'<text x="{x_right}" y="{y + 72}">Outlined: {"cell" if len(highlight) == 1 else "cells"} '
                   f"{marks} (each in its colour, in both panels)</text>")
    notes = (f"Exact data: {name} (SHA-256 {digest[:16]}).",
             "Coordinates are rounded to 0.001 px when drawn; the certificate, not this figure, is the "
             "mathematical object.",
             ("Labels: cell indices and point indices k of a_k ↦ b_k (0-based). " if labels else
              "Labels: point indices k of a_k ↦ b_k (0-based). ") + "No layer order is drawn or implied.")
    for n, note in enumerate(notes):
        svg.append(f'<text x="{LEFT}" y="{height - 40 + 14 * n}" font-size="10" fill="#404040">{escape(note)}</text>')
    svg.append("</svg>")
    return "\n".join(svg) + "\n", {"cells": m, "faces": nfaces, "creases": count["crease"],
                                   "auxiliary": count["aux"], "boundary": count["edge"]}


# --- command line ------------------------------------------------------------------------------

def run_check(exe, path):
    """brehm_cpp check on the very bytes that will be drawn. Returns (exit code, report or message)."""
    exe = str(Path(exe).resolve()) if Path(exe).is_file() else exe  # Windows cannot start "dir/exe" relative paths
    try:
        r = subprocess.run([exe, "check", str(path)], capture_output=True, text=True)
    except OSError as e:
        return 2, f"cannot run {exe}: {e}"
    if r.returncode == 0:
        return 0, json.loads(r.stdout)
    if r.returncode == 1:
        try:
            why = "; ".join(json.loads(r.stdout).get("failures", []))
        except ValueError:
            why = r.stderr.strip()
        return 1, "certificate rejected by brehm_cpp check: " + why[:2000]
    return {3: 3, 5: 5}.get(r.returncode, 4), f"brehm_cpp check exited {r.returncode}: {r.stderr.strip()}"


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 epilog="Exit codes: 0 written, 1 rejected, 2 usage, 3 unreadable input, "
                                        "4 internal error, 5 out of memory.")
    ap.add_argument("cert", help="certificate (brehm-rational-v1)")
    ap.add_argument("out", help="figure to write (.svg)")
    how = ap.add_mutually_exclusive_group()
    how.add_argument("--exe", help="brehm_cpp executable used to check the certificate (default: on PATH)")
    how.add_argument("--unverified", action="store_true", help="draw without checking; the figure says so")
    lab = ap.add_mutually_exclusive_group()
    lab.add_argument("--labels", dest="labels", action="store_true", default=None,
                     help=f"label every cell with its index (default up to {AUTO_LABELS} cells)")
    lab.add_argument("--no-labels", dest="labels", action="store_false")
    ap.add_argument("--highlight", default="", metavar="I,J,...", help="outline these cells in both panels")
    ap.add_argument("--zoom", default="", metavar="I,J,...",
                    help="fit both panels to these cells (the others are still drawn, clipped)")
    ap.add_argument("--no-points", action="store_true", help="omit the sources and targets")
    args = ap.parse_args(argv)
    cert_path, out = Path(args.cert), Path(args.out)
    if out.suffix.lower() != ".svg":
        print("the output file must end in .svg", file=sys.stderr)
        return 2
    if out.exists() and cert_path.exists() and os.path.samefile(out, cert_path):
        print("the output file is the certificate itself", file=sys.stderr)
        return 2
    try:
        highlight, zoom = ([int(k) for k in arg.split(",") if k.strip()] for arg in (args.highlight, args.zoom))
    except ValueError:
        print("--highlight and --zoom take cell indices, e.g. 3,7", file=sys.stderr)
        return 2
    exe = None
    if not args.unverified:
        exe = args.exe or shutil.which("brehm_cpp")
        if exe is None:
            print("brehm_cpp not found on PATH: pass --exe PATH, or --unverified to draw without checking",
                  file=sys.stderr)
            return 2
    t0 = time.perf_counter()
    try:
        data = cert_path.read_bytes()
        digest = hashlib.sha256(data).hexdigest()
        with tempfile.TemporaryDirectory() as tmp:
            copy = Path(tmp) / "cert.json"
            copy.write_bytes(data)
            report = None
            if exe:
                code, report = run_check(exe, copy)
                if code:
                    print(report, file=sys.stderr)
                    return code
            raw = load(copy)
        del data
        domain, sources, targets, cells = read(raw)
        if any(not 0 <= k < len(cells) for k in highlight + zoom):
            print(f"--highlight, --zoom: cell indices run from 0 to {len(cells) - 1}", file=sys.stderr)
            return 2
        intervals = structure(domain, cells)
        name = cert_path.name
        labels = len(cells) <= AUTO_LABELS if args.labels is None else args.labels
        text, summary = draw(name, digest, domain, sources, targets, cells, intervals, raw.get("construction"),
                             report is not None, labels, highlight, not args.no_points, zoom)
        if report is not None:  # the checker's independent count of the same structure
            mine = {"pieces": len(cells), "interior_intervals": summary["creases"] + summary["auxiliary"],
                    "boundary_intervals": summary["boundary"], "isometric_regions": summary["faces"],
                    "distinct_maps": len({c[1:] for c in cells})}
            theirs = {k: report.get(k) for k in mine}
            if mine != theirs:
                print(f"internal error: edge structure {mine} differs from brehm_cpp check {theirs}", file=sys.stderr)
                return 4
    except Rejected as e:
        print(f"certificate rejected: {e}", file=sys.stderr)
        return 1
    except (OSError, ValueError, KeyError, TypeError, AttributeError, IndexError, RecursionError) as e:
        print(f"invalid input: {e!r}", file=sys.stderr)
        return 3
    except MemoryError:
        print("out of memory", file=sys.stderr)
        return 5
    except Exception as e:  # noqa: BLE001 - any other failure is ours, not a verdict on the input
        print(f"internal error: {e!r}", file=sys.stderr)
        return 4
    try:
        part = out.with_name(out.name + ".part")
        part.write_text(text, encoding="utf-8", newline="\n")
        os.replace(part, out)
    except OSError as e:
        print(f"cannot write {out}: {e}", file=sys.stderr)
        return 4
    summary.update(svg=out.name, verified=report is not None, bytes=len(text.encode("utf-8")),
                   seconds=round(time.perf_counter() - t0, 3))
    print(json.dumps(summary))
    return 0


if __name__ == "__main__":
    sys.set_int_max_str_digits(0)  # long exact numbers are legitimate
    sys.exit(main())
