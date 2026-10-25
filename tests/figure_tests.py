"""Independent tests of tools/brehm_figure.py: a figure must show exactly its certificate.

    python tests/figure_tests.py [path/to/brehm_cpp] [--random N] [--tool path/to/brehm_figure.py]

Every example (and N random instances from tests/fuzz.py) is built with brehm_cpp and drawn. The SVG
is then parsed, and every drawn element is recomputed here from the certificate with this file's own
exact arithmetic (no code is shared with the renderer):

* each cell is drawn once in each panel, its polygon on the left and its image A v + b on the right,
  every vertex the exact point rounded to 0.001 px; the left tint is the sign of det A, the hatching
  on the right marks exactly the cells with det A = -1, and all image cells share one fill;
* the drawing transform is one uniform scale with the y axis flipped, the same in both panels, the
  boxes are centred and the larger one fills its panel;
* edges: every pair of cell edges is intersected here (not the renderer's grouping by lines); a common
  piece with equal maps is auxiliary, with different maps a crease (both maps must agree on it), a
  piece on the paper's boundary is boundary. The segments drawn for each kind in each panel must be
  exactly these, and their numbers must equal brehm_cpp check's interval and region counts;
* sources, targets and cell labels at their exact positions.

Then: the figures of the two stored certificates must be byte for byte the recorded ones (a visual
regression test; the output is exact and deterministic); identical geometry after scaling by 10^30
and 10^-30 and after translations; rotations and reflections; cell order and vertex rotation;
byte-identical output across runs and hash seeds; the certificate is never modified; rejected or malformed input gives no figure and the documented exit
code; the renderer's cross-check against brehm_cpp check catches a dropped edge. Exit code 1 on any
failure.
"""
import copy
import hashlib
import json
import os
import random
import re
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from collections import Counter
from fractions import Fraction as F
from pathlib import Path

CPP = Path(__file__).resolve().parents[1]
TOOL = CPP / "tools" / "brehm_figure.py"
SVG = "{http://www.w3.org/2000/svg}"
FRONT, BACK, INK = "#f4ecd2", "#c9dbee", "#1b2a41"            # the documented colours of the figure
KIND = {"#808080": "aux", "#000000": "edge", "#c0392b": "crease"}
PAIRWISE_MAX = 400                                           # pairwise edge test up to this many cells
GOLDEN = {  # SHA-256 of the figures of the stored certificates: a byte-level visual regression test.
            # The output is exact and deterministic, so any change to a drawing shows here; update
            # these only for an intended change, after checking the new figures.
    "regression_line_stored": "5afa592d87f826c37d4b9d405e3af96a77ab6cd8ec45acf880d7093d5d7ccd99",
    "regression_fold_stored": "043e5d3e65eba86e2f956d0c6c5e827d400f322c35842417b31fe2840353cdc9"}
failures = []


def expect(cond, what):
    if not cond:
        failures.append(what)
        print("FAIL", what, flush=True)
    return cond


# --- exact geometry, written independently of the renderer --------------------------------------

def load_cert(path):
    raw = json.loads(Path(path).read_text(encoding="utf-8"))
    P = lambda v: (F(v[0]), F(v[1]))  # noqa: E731
    cells = [([P(v) for v in c["polygon"]], (P(c["A"][0]), P(c["A"][1])), P(c["b"])) for c in raw["pieces"]]
    return [P(v) for v in raw["domain"]], [P(v) for v in raw["sources"]], [P(v) for v in raw["targets"]], cells


def img(cell, p):
    (a, b), (c, d) = cell[1]
    return (a * p[0] + b * p[1] + cell[2][0], c * p[0] + d * p[1] + cell[2][1])


def det(cell):
    (a, b), (c, d) = cell[1]
    return a * d - b * c


def cross(o, a, b):
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def edges(poly):
    return [(poly[i], poly[(i + 1) % len(poly)]) for i in range(len(poly))]


def common(e, f):
    """The common piece of positive length of two segments, or None."""
    (p, q), (r, t) = e, f
    if cross(p, q, r) or cross(p, q, t):
        return None
    d = (q[0] - p[0], q[1] - p[1])
    dd = d[0] ** 2 + d[1] ** 2
    par = lambda x: ((x[0] - p[0]) * d[0] + (x[1] - p[1]) * d[1]) / dd  # noqa: E731
    lo, hi = max(F(0), min(par(r), par(t))), min(F(1), max(par(r), par(t)))
    return ((p[0] + lo * d[0], p[1] + lo * d[1]), (p[0] + hi * d[0], p[1] + hi * d[1])) if lo < hi else None


def bbox(points):
    return (min(p[0] for p in points), min(p[1] for p in points), max(p[0] for p in points), max(p[1] for p in points))


def pieces(domain, cells):
    """[(kind, i, j, segment)] by intersecting every pair of edges of cells with overlapping boxes."""
    boxes, E, out = [bbox(c[0]) for c in cells], [edges(c[0]) for c in cells], []
    for i in range(len(cells)):
        for j in range(i + 1, len(cells)):
            (a, b, c, d), (e, f, g, h) = boxes[i], boxes[j]
            if e > c or a > g or f > d or b > h:
                continue
            for x in E[i]:
                for y in E[j]:
                    seg = common(x, y)
                    if seg:
                        expect(all(img(cells[i], p) == img(cells[j], p) for p in seg), f"cells {i}, {j} disagree")
                        out.append(("aux" if cells[i][1:] == cells[j][1:] else "crease", i, j, seg))
        for x in E[i]:
            for y in edges(domain):
                seg = common(x, y)
                if seg:
                    out.append(("edge", i, None, seg))
    return out


def regions(n, found):
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            i = parent[i]
        return i

    for kind, i, j, _ in found:
        if kind == "aux":
            parent[find(i)] = find(j)
    return sum(find(i) == i for i in range(n))


def r3(q):
    """The nearest multiple of 1/1000, ties to even."""
    return F(round(q * 1000), 1000)


# --- reading a figure ---------------------------------------------------------------------------

def pts(text):
    return [tuple(F(c) for c in xy.split(",")) for xy in text.split()]


def parse(path):
    root = ET.parse(path).getroot()
    fig = {"desc": json.loads(root.find(SVG + "desc").text)}
    for cp in root.iter(SVG + "clipPath"):
        r = cp.find(SVG + "rect")
        fig[cp.get("id")] = (tuple(F(r.get(k)) for k in ("x", "y", "width", "height")) if r is not None
                             else pts(cp.find(SVG + "polygon").get("points")))
    for name in ("paper", "image"):
        g = root.find(f"{SVG}g[@id='{name}']")
        panel = {"tint": {}, "fill": [], "hatch": [], "segments": Counter(), "outlines": [], "circles": [],
                 "point_labels": [], "cell_labels": [], "fill_groups": []}
        for el in g:
            tag = el.tag[len(SVG):]
            if tag == "g":
                fill = el.get("fill")
                polys = list(el.iter(SVG + "polygon"))
                if fill in (FRONT, BACK):
                    for p in polys:
                        panel["tint"][p.find(SVG + "title").text] = (fill, pts(p.get("points")))
                elif fill == "url(#reversed)":
                    panel["hatch"] += [pts(p.get("points")) for p in polys]
                else:
                    panel["fill_groups"].append((fill, el.get("fill-opacity"),
                                                 [p.get("fill") or p.get("fill-opacity") for p in polys]))
                    panel["fill"] += [(p.find(SVG + "title").text, pts(p.get("points"))) for p in polys]
            elif tag == "path":
                kind = KIND[el.get("stroke")]
                for m in re.finditer(r"M(\S+) (\S+)L(\S+) (\S+?)(?=M|$)", el.get("d")):
                    a, b = (F(m.group(1)), F(m.group(2))), (F(m.group(3)), F(m.group(4)))
                    panel["segments"][(kind, min(a, b), max(a, b))] += 1
            elif tag == "polygon":
                panel["outlines"].append((el.get("stroke"), el.get("clip-path"), pts(el.get("points"))))
            elif tag == "circle":
                panel["circles"].append((F(el.get("cx")), F(el.get("cy"))))
            elif tag == "text":
                where = (F(el.get("x")), F(el.get("y")))
                (panel["cell_labels"] if el.get("font-size") == "9" else panel["point_labels"]).append((where, el.text))
        fig[name] = panel
    return fig


def indices(label):
    out = []
    for part in label.split(","):
        a, _, b = part.partition("–")
        out += list(range(int(a), int(b) + 1)) if b else [int(a)]
    return out


# --- the checks ---------------------------------------------------------------------------------

def verify(cert_path, svg_path, report, zoom=(), labels=None, points=True, tag=""):
    """Every check of the module docstring for one figure. Returns the number of failures added."""
    before = len(failures)
    domain, sources, targets, cells = load_cert(cert_path)
    fig, m = parse(svg_path), len(cells)
    desc = fig["desc"]
    s = F(desc["scale_px_per_unit"])
    views = {"paper": tuple(map(F, desc["paper_transform"])), "image": tuple(map(F, desc["image_transform"]))}
    show = lambda name, p: (r3(views[name][0] + s * p[0]), r3(views[name][1] - s * p[1]))  # noqa: E731
    images = [[img(c, p) for p in c[0]] for c in cells]
    expect(desc["cells"] == m, f"{tag}: desc says {desc['cells']} cells, certificate has {m}")
    # transform: uniform scale, y flipped, boxes centred, the larger box fills its panel
    if not zoom:
        boxes = {"paper": bbox(domain + (sources if points else [])),
                 "image": bbox([p for poly in images for p in poly] + (targets if points else []))}
        span = max(max(b[2] - b[0], b[3] - b[1]) for b in boxes.values())
        for name, clip in (("paper", fig["paper-panel"]), ("image", fig["image-panel"])):
            x0, y0, w, h = clip
            xmin, ymin, xmax, ymax = boxes[name]
            expect(s == w / span and w == h, f"{tag}: scale {s} is not panel/span = {w / span}")
            expect(views[name][0] + s * xmin - x0 == x0 + w - (views[name][0] + s * xmax), f"{tag}: {name} not centred in x")
            expect(views[name][1] - s * ymax - y0 == y0 + h - (views[name][1] - s * ymin), f"{tag}: {name} not centred in y")
    # cells: once per panel, exact vertices, tint = orientation, one common fill, hatching = reversed
    left, right = fig["paper"], fig["image"]
    expect(len(left["tint"]) == m and len(right["fill"]) == m, f"{tag}: {len(left['tint'])} / {len(right['fill'])} cells drawn, not {m}")
    expect(len(right["fill_groups"]) == 1 and right["fill_groups"][0][0] == INK
           and all(x is None for x in right["fill_groups"][0][2]), f"{tag}: image cells do not share one fill")
    for k, c in enumerate(cells):
        title = f"cell {k}" + ("" if det(c) > 0 else " (reversed)")
        tint, drawn = left["tint"].get(title, (None, None))
        expect(drawn == [show("paper", p) for p in c[0]], f"{tag}: cell {k} drawn wrongly on the paper")
        expect(tint == (FRONT if det(c) > 0 else BACK), f"{tag}: cell {k} has the wrong tint")
        expect(det(c) in (1, -1), f"{tag}: cell {k} has det {det(c)}")
    for k, (title, drawn) in enumerate(right["fill"]):
        expect(title.split()[1] == str(k) and drawn == [show("image", p) for p in images[k]],
               f"{tag}: image of cell {k} drawn wrongly")
    expect(Counter(map(tuple, right["hatch"])) == Counter(tuple(show("image", p) for p in images[k])
                                                          for k in range(m) if det(cells[k]) < 0),
           f"{tag}: hatching does not mark exactly the reversed cells")
    # edges: pairwise, independently of the renderer; then the counts of brehm_cpp check
    if m <= PAIRWISE_MAX:
        found = pieces(domain, cells)
        for name, panel in (("paper", left), ("image", right)):
            want = Counter()
            for kind, i, _, seg in found:
                a, b = (show(name, img(cells[i], p) if name == "image" else p) for p in seg)
                want[(kind, min(a, b), max(a, b))] += 1
            expect(panel["segments"] == want, f"{tag}: {name} edges differ: {sum((panel['segments'] - want).values())} "
                                              f"extra, {sum((want - panel['segments']).values())} missing")
        if report:
            n = Counter(k for k, *_ in found)
            mine = (n["crease"] + n["aux"], n["edge"], regions(m, found))
            theirs = (report["interior_intervals"], report["boundary_intervals"], report["isometric_regions"])
            expect(mine == theirs, f"{tag}: pieces/regions {mine} differ from brehm_cpp check {theirs}")
    # points and labels
    if points:
        for name, panel, data in (("paper", left, sources), ("image", right, targets)):
            groups = {}
            for k, p in enumerate(data):
                groups.setdefault(p, []).append(k)
            want = sorted((show(name, p), ks) for p, ks in groups.items())
            expect(sorted(panel["circles"]) == [w for w, _ in want], f"{tag}: {name} points misplaced")
            expect(sorted((w, indices(t)) for w, t in panel["point_labels"]) == want, f"{tag}: {name} point labels wrong")
    if labels if labels is not None else m <= 60:
        for name, panel in (("paper", left), ("image", right)):
            want = {}
            for k, c in enumerate(cells):
                mid = (sum(p[0] for p in c[0]) / len(c[0]), sum(p[1] for p in c[0]) / len(c[0]))
                want.setdefault(show(name, img(c, mid) if name == "image" else mid), []).append(k)
            got = {w: indices(t) for w, t in panel["cell_labels"]}
            expect(got == want, f"{tag}: {name} cell labels wrong")
    return len(failures) - before


# --- running the renderer -----------------------------------------------------------------------

def draw(cert, out, *opts, env=None):
    r = subprocess.run([sys.executable, str(TOOL), str(cert), str(out), *map(str, opts)],
                       capture_output=True, text=True, env=env)
    return r.returncode, r.stdout, r.stderr


def check_report(cert):
    r = subprocess.run([EXE, "check", str(cert)], capture_output=True, text=True)
    return json.loads(r.stdout) if r.returncode == 0 else None


def geometry(svg_path):
    """Everything drawn inside the two panels, in document order (no text)."""
    root = ET.parse(svg_path).getroot()
    out = []
    for name in ("paper", "image"):
        for el in root.find(f"{SVG}g[@id='{name}']").iter():
            if el.tag in (SVG + "polygon", SVG + "path", SVG + "circle"):
                out.append((el.tag, el.get("points") or el.get("d") or (el.get("cx"), el.get("cy"))))
    return out


def transformed(raw, src=None, tgt=None):
    """The certificate of x -> N(F(M^-1 x)) on M(D), for affine maps M (source side) and N (target
    side), each given as (2x2 matrix, vector); a scaling by lam is ((lam, 0), (0, lam)), 0."""
    c = copy.deepcopy(raw)
    P = lambda v: (F(v[0]), F(v[1]))  # noqa: E731
    S = lambda p: [str(p[0]), str(p[1])]  # noqa: E731
    mv = lambda M, p: (M[0][0] * p[0] + M[0][1] * p[1], M[1][0] * p[0] + M[1][1] * p[1])  # noqa: E731
    mm = lambda X, Y: tuple(tuple(sum(X[i][k] * Y[k][j] for k in range(2)) for j in range(2)) for i in range(2))  # noqa: E731
    I = ((F(1), F(0)), (F(0), F(1)))
    (M, t), (N, u) = src or (I, (0, 0)), tgt or (I, (0, 0))
    dM = M[0][0] * M[1][1] - M[0][1] * M[1][0]
    Minv = ((M[1][1] / dM, -M[0][1] / dM), (-M[1][0] / dM, M[0][0] / dM))
    fwd = lambda p: (mv(M, p)[0] + t[0], mv(M, p)[1] + t[1])  # noqa: E731
    c["domain"] = [S(fwd(P(v))) for v in c["domain"]]
    c["sources"] = [S(fwd(P(v))) for v in c["sources"]]
    c["targets"] = [S((mv(N, P(v))[0] + u[0], mv(N, P(v))[1] + u[1])) for v in c["targets"]]
    for piece in c["pieces"]:
        A = (P(piece["A"][0]), P(piece["A"][1]))
        b = P(piece["b"])
        A2 = mm(N, mm(A, Minv))                    # F'(y) = N A M^-1 (y - t) + N b + u
        At = mv(A2, t)
        Nb = mv(N, b)
        piece["polygon"] = [S(fwd(P(v))) for v in piece["polygon"]]
        piece["A"] = [S(A2[0]), S(A2[1])]
        piece["b"] = S((Nb[0] + u[0] - At[0], Nb[1] + u[1] - At[1]))
    return c


def instances(tmp, n_random):
    """(name, certificate path) for every buildable example and n_random fuzz instances."""
    out = []
    for inst in sorted((CPP / "examples").glob("*.json")):
        if inst.name.endswith("_cert.json") or json.loads(inst.read_text())["expect"] != "ok":
            continue
        cert = Path(tmp) / f"{inst.stem}_cert.json"
        r = subprocess.run([EXE, "construct", str(inst), str(cert)], capture_output=True, text=True)
        if expect(r.returncode == 0, f"{inst.stem}: construct exit {r.returncode}"):
            out.append((inst.stem, cert))
    (Path(tmp) / "stored").mkdir()
    for name in ("regression_line", "regression_fold"):  # with LF endings, as released: a figure names its file's SHA-256
        copy_ = Path(tmp) / "stored" / f"{name}_cert.json"
        copy_.write_bytes((CPP / "examples" / f"{name}_cert.json").read_bytes().replace(b"\r\n", b"\n"))
        out.append((name + "_stored", copy_))
    sys.path.insert(0, str(CPP / "tests"))
    import fuzz  # noqa: E402
    rng = random.Random(11)
    kinds = ["contraction", "rotation_contraction", "projection", "fold", "collapse", "rigid", "identity"]
    for it in range(n_random):
        dom = fuzz.paper(rng)
        n = rng.randint(1, 6)
        src = fuzz.sources(rng, dom, n)
        tgt = fuzz.targets(rng, src, rng.choice(kinds))
        order = list(range(n))
        rng.shuffle(order)
        c, s = rng.choice(fuzz.PYTH)
        g = [[c, -s], [s, c]] if rng.random() < 0.5 else [[c, s], [s, -c]]
        st = lambda p: [str(p[0]), str(p[1])]  # noqa: E731
        inst = {"domain": [st(p) for p in dom], "sources": [st(p) for p in src], "targets": [st(p) for p in tgt],
                "order": order, "g": [[str(x) for x in r] for r in g], "L": str(fuzz.admissible_L(rng, dom, src, src[order[0]]))}
        ip, cp = Path(tmp) / f"random{it}.json", Path(tmp) / f"random{it}_cert.json"
        ip.write_text(json.dumps(inst))
        r = subprocess.run([EXE, "construct", str(ip), str(cp)], capture_output=True, text=True)
        if expect(r.returncode == 0, f"random {it}: construct exit {r.returncode}"):
            out.append((f"random{it}", cp))
    return out


def main():
    global EXE, TOOL
    args = sys.argv[1:]
    n_random = 0
    if "--random" in args:
        i = args.index("--random")
        n_random = int(args[i + 1])
        del args[i:i + 2]
    if "--tool" in args:  # test another copy of the renderer
        i = args.index("--tool")
        TOOL = Path(args[i + 1]).resolve()
        del args[i:i + 2]
    EXE = args[0] if args else str(next(p for p in (CPP / "build" / "brehm_cpp.exe", CPP / "build" / "brehm_cpp")
                                        if p.exists()))
    EXE = str(Path(EXE).resolve()) if Path(EXE).is_file() else EXE  # Windows cannot start "dir/exe" relative paths
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        certs = instances(tmp, n_random)
        for name, cert in certs:  # 1. every figure against its certificate
            out = tmp / f"{name}.svg"
            code, stdout, err = draw(cert, out, "--exe", EXE)
            if not expect(code == 0, f"{name}: exit {code} {err.strip()}"):
                continue
            n = verify(cert, out, check_report(cert), tag=name)
            if name in GOLDEN:
                n += not expect(hashlib.sha256(out.read_bytes()).hexdigest() == GOLDEN[name],
                                f"{name}: the figure differs from the recorded one (visual regression)")
            print(f"{name:40s} {json.loads(stdout)['cells']:5d} cells  {'ok' if not n else f'{n} FAILURES'}", flush=True)

        base = tmp / "regression_line_cert.json"
        raw = json.loads(base.read_text())
        draw(base, tmp / "base.svg", "--exe", EXE)
        ref = geometry(tmp / "base.svg")
        # 2. scaling and translation leave the drawing unchanged; rotation and reflection pass every check
        I = ((F(1), F(0)), (F(0), F(1)))
        cases = [("scaled by 10^30", (((F(10**30), 0), (0, F(10**30))), (0, 0)), (((F(10**30), 0), (0, F(10**30))), (0, 0)), True),
                 ("scaled by 10^-30", (((F(1, 10**30), 0), (0, F(1, 10**30))), (0, 0)), (((F(1, 10**30), 0), (0, F(1, 10**30))), (0, 0)), True),
                 ("translated", (I, (F(10**12, 7), F(-3, 11))), (I, (F(-5, 3), F(10**9))), True),
                 ("paper rotated by (3/5, 4/5)", (((F(3, 5), F(-4, 5)), (F(4, 5), F(3, 5))), (0, 0)), None, False),
                 ("image reflected", None, (((F(1), F(0)), (F(0), F(-1))), (F(1, 3), 0)), False),
                 ("both rotated, far away", (((F(5, 13), F(-12, 13)), (F(12, 13), F(5, 13))), (F(10**20), 0)),
                  (((F(-8, 17), F(15, 17)), (F(15, 17), F(8, 17))), (0, F(-10**20))), False)]
        for what, src, tgt, same in cases:
            path, out = tmp / "moved.json", tmp / "moved.svg"
            path.write_text(json.dumps(transformed(raw, src, tgt)))
            report = check_report(path)
            if not expect(report is not None, f"{what}: brehm_cpp check rejects the moved certificate"):
                continue
            code, _, err = draw(path, out, "--exe", EXE)
            if expect(code == 0, f"{what}: exit {code} {err.strip()}"):
                n = verify(path, out, report, tag=what)
                if same:
                    expect(geometry(out) == ref, f"{what}: drawing differs from the original")
                print(f"{what:40s} {'ok' if not n else f'{n} FAILURES'}", flush=True)
        # 3. cell order and vertex rotation: the same drawing up to the order of elements
        shuffled = copy.deepcopy(raw)
        shuffled["pieces"].reverse()
        for p in shuffled["pieces"]:
            p["polygon"].append(p["polygon"].pop(0))
        (tmp / "shuffled.json").write_text(json.dumps(shuffled))
        draw(tmp / "shuffled.json", tmp / "shuffled.svg", "--exe", EXE)
        verify(tmp / "shuffled.json", tmp / "shuffled.svg", check_report(tmp / "shuffled.json"), tag="shuffled")
        norm = lambda poly: min(tuple(poly[i:] + poly[:i]) for i in range(len(poly)))  # noqa: E731
        for name in ("paper", "image"):
            a, b = parse(tmp / "base.svg")[name], parse(tmp / "shuffled.svg")[name]
            expect(a["segments"] == b["segments"], f"shuffled: {name} edges differ")
            cells_a = sorted(norm(p) for _, p in (a["tint"].values() if name == "paper" else a["fill"]))
            cells_b = sorted(norm(p) for _, p in (b["tint"].values() if name == "paper" else b["fill"]))
            expect(cells_a == cells_b, f"shuffled: {name} cells differ")
        print("cell order and vertex rotation                 done", flush=True)
        # 4. determinism, and the certificate is never written
        digest = hashlib.sha256(base.read_bytes()).hexdigest()
        outs = []
        for seed in ("0", "1", "2"):
            env = dict(os.environ, PYTHONHASHSEED=seed)
            draw(base, tmp / f"det{seed}.svg", "--exe", EXE, "--labels", "--highlight", "3,1", env=env)
            outs.append((tmp / f"det{seed}.svg").read_bytes())
        expect(len(set(outs)) == 1, "output differs between runs or hash seeds")
        expect(hashlib.sha256(base.read_bytes()).hexdigest() == digest, "the certificate was modified")
        # 5. options: zoom, highlight, no points, no labels, unverified
        code, _, _ = draw(base, tmp / "zoom.svg", "--exe", EXE, "--zoom", "14,8", "--highlight", "14,8")
        if expect(code == 0, "zoom: exit code"):
            verify(base, tmp / "zoom.svg", check_report(base), zoom=(14, 8), tag="zoom")
            fig = parse(tmp / "zoom.svg")
            polys = {"paper": {int(t.split()[1]): p for t, (_, p) in fig["paper"]["tint"].items()},
                     "image": {int(t.split()[1]): p for t, p in fig["image"]["fill"]}}
            for name in ("paper", "image"):  # a hairline, then a band clipped to the cell itself
                got = fig[name]["outlines"]
                expect([p for _, _, p in got] == [polys[name][14], polys[name][8]] * 2,
                       f"highlight outlines are not the cells' polygons ({name})")
                expect([clip for _, clip, _ in got[:2]] == [None, None] and
                       all(clip and fig[clip[5:-1]] == p for _, clip, p in got[2:]),
                       f"highlight bands are not clipped to their own cells ({name})")
            expect([c for c, _, _ in fig["paper"]["outlines"]] == [c for c, _, _ in fig["image"]["outlines"]],
                   "a highlighted cell has different colours in the two panels")
        code, _, _ = draw(base, tmp / "plain.svg", "--exe", EXE, "--no-points", "--no-labels")
        if expect(code == 0, "no points: exit code"):
            verify(base, tmp / "plain.svg", check_report(base), labels=False, points=False, tag="no points")
        code, _, _ = draw(base, tmp / "unverified.svg", "--unverified")
        expect(code == 0 and "NOT VERIFIED" in (tmp / "unverified.svg").read_text(encoding="utf-8"), "unverified: not marked")
        # 6. failures: no figure, documented exit code
        bad = copy.deepcopy(raw)
        bad["pieces"][5]["b"][1] = str(F(bad["pieces"][5]["b"][1]) + F(1, 1000))
        dup = copy.deepcopy(raw)
        dup["pieces"].append(copy.deepcopy(dup["pieces"][2]))
        (tmp / "bad.json").write_text(json.dumps(bad))
        (tmp / "dup.json").write_text(json.dumps(dup))
        (tmp / "broken.json").write_text('{"format": "brehm-rational-v1", "domain": [')
        (tmp / "twice.json").write_text(raw_text := base.read_text().replace('"format"', '"format":"x","format"', 1))
        (tmp / "deep.json").write_text('{"format": "brehm-rational-v1", "x": ' + "[" * 70 + "]" * 70 + "}")
        (tmp / "fake.svg").write_text(base.read_text())
        runs = [("rejected by brehm_cpp check", ("bad.json", "o1.svg", "--exe", EXE), 1),
                ("cells overlap (unverified)", ("dup.json", "o2.svg", "--unverified"), 1),
                ("truncated JSON", ("broken.json", "o3.svg", "--unverified"), 3),
                ("duplicate key", ("twice.json", "o4.svg", "--exe", EXE), 3),
                ("duplicate key (unverified)", ("twice.json", "o5.svg", "--unverified"), 3),
                ("nesting deeper than 64", ("deep.json", "o6.svg", "--unverified"), 3),
                ("output is not .svg", ("regression_line_cert.json", "o7.png", "--exe", EXE), 2),
                ("output is the certificate", ("fake.svg", "fake.svg", "--unverified"), 2),
                ("no such executable", ("regression_line_cert.json", "o8.svg", "--exe", str(tmp / "nope")), 2),
                ("cell index out of range", ("regression_line_cert.json", "o9.svg", "--exe", EXE, "--highlight", str(len(raw["pieces"]))), 2),
                ("missing file", ("nothing.json", "o10.svg", "--unverified"), 3)]
        for what, (cert, out, *opts), want in runs:
            code, _, err = draw(tmp / cert, tmp / out, *opts)
            made = (tmp / out).exists() and cert != out
            expect(code == want and not made and not (tmp / (out + ".part")).exists(),
                   f"{what}: exit {code} (expected {want}), figure written: {made}; {err.strip()[:200]}")
        expect((tmp / "fake.svg").read_text() == base.read_text(), "the certificate named .svg was overwritten")
        # 7. the renderer's own cross-check: a dropped edge piece must stop it (exit 4, no figure)
        sys.path.insert(0, str(TOOL.parent))
        import brehm_figure  # noqa: E402
        real = brehm_figure.structure
        brehm_figure.structure = lambda d, c: real(d, c)[:-1]
        code = brehm_figure.main([str(base), str(tmp / "dropped.svg"), "--exe", EXE])
        brehm_figure.structure = real
        expect(code == 4 and not (tmp / "dropped.svg").exists(), f"dropped edge piece: exit {code}")
        print("options, failures and the cross-check          done", flush=True)
    print("ALL PASSED" if not failures else f"{len(failures)} FAILURES")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.set_int_max_str_digits(0)
    sys.exit(main())
