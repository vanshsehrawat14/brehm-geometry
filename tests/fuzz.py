"""Seeded random testing of brehm_cpp with all construction options, aimed at equality branches.

    python tests/fuzz.py N SEED [path/to/brehm_cpp]

Each instance: a random convex rational paper (hull of grid points), 1-6 sources drawn from its
vertices, edge points, grid points and repeats (so cut lines through vertices, sources on the
boundary and duplicates occur), short data of one of seven kinds (contraction, rotation-contraction,
projection, fold, collapse, rigid motion: every pair tight, identity), a random insertion order, a
random exact rotation or reflection g, and a random admissible L from just above Brehm's bound to
100 times it. Every output must pass brehm_cpp check, tools/ref_check.py --instance and agree with
tools/brehm_eval.py at every cell vertex. Prints a summary line; exit code 1 on any failure.
"""
import json
import random
import subprocess
import sys
import tempfile
from fractions import Fraction as F
from pathlib import Path

CPP = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(CPP / "tools"))
import brehm_eval  # noqa: E402
import ref_check  # noqa: E402

PYTH = [(F(1), F(0)), (F(3, 5), F(4, 5)), (F(5, 13), F(12, 13)), (F(8, 17), F(15, 17)), (F(0), F(1)), (F(-1), F(0))]


def hull(pts):
    pts = sorted(set(pts))
    if len(pts) < 3:
        return pts

    def half(seq):
        out = []
        for p in seq:
            while len(out) >= 2 and ref_check.turn(out[-2], out[-1], p) <= 0:
                out.pop()
            out.append(p)
        return out[:-1]
    return half(pts) + half(pts[::-1])


def paper(rng):
    while True:
        h = hull([(F(rng.randint(-8, 8), 8), F(rng.randint(-8, 8), 8)) for _ in range(rng.randint(3, 9))])
        if len(h) >= 3:
            return h


def sources(rng, dom, n):
    out = []
    while len(out) < n:
        r = rng.random()
        if r < 0.2:
            p = rng.choice(dom)
        elif r < 0.4:
            u, v = rng.choice([(dom[i], dom[(i + 1) % len(dom)]) for i in range(len(dom))])
            t = F(rng.randint(0, 4), 4)
            p = (u[0] + t * (v[0] - u[0]), u[1] + t * (v[1] - u[1]))
        elif r < 0.5 and out:
            p = rng.choice(out)
        else:
            w = [F(rng.randint(0, 4)) for _ in dom]
            s = sum(w) or F(1)
            p = (sum(wi * v[0] for wi, v in zip(w, dom)) / s, sum(wi * v[1] for wi, v in zip(w, dom)) / s)
        if ref_check.inside(dom, p):
            out.append(p)
    return out


def targets(rng, pts, kind):
    c, s = rng.choice(PYTH)
    shift = (F(rng.randint(-4, 4), 4), F(rng.randint(-4, 4), 4))
    rot = lambda p: (c * p[0] - s * p[1], s * p[0] + c * p[1])  # noqa: E731
    lam = rng.choice([F(1, 2), F(1, 3), F(2, 3)])
    if kind == "contraction":
        return [(lam * x, lam * y) for x, y in pts]
    if kind == "rotation_contraction":
        return [rot((lam * x, lam * y)) for x, y in pts]
    if kind == "projection":
        return [(x, 0) for x, _ in pts]
    if kind == "fold":
        return [(abs(x - shift[0]), y) for x, y in pts]
    if kind == "collapse":
        return [shift] * len(pts)
    if kind == "rigid":
        return [(rot(p)[0] + shift[0], rot(p)[1] + shift[1]) for p in pts]
    return list(pts)


def admissible_L(rng, dom, pts, a1):
    need = 2 * max(((p[0] - a1[0]) ** 2 + (p[1] - a1[1]) ** 2 for p in pts), default=F(0))
    reach = max(max(abs(p[0] - a1[0]), abs(p[1] - a1[1])) for p in dom) / 3
    L = max(F(int(float(need) ** 0.5 * 1000) + 1, 1000), reach, F(1, 1000))
    while L * L < need:
        L += F(1, 1000)
    return L * rng.choice([F(1), F(101, 100), F(3, 2), F(4), F(100)])


def main():
    N, seed = int(sys.argv[1]), int(sys.argv[2])
    exe = sys.argv[3] if len(sys.argv) > 3 else str(next(p for p in (CPP / "build" / "brehm_cpp.exe",
                                                                       CPP / "build" / "brehm_cpp") if p.exists()))
    exe = str(Path(exe).resolve()) if Path(exe).is_file() else exe  # Windows cannot start "dir/exe" relative paths
    rng = random.Random(seed)
    kinds = ["contraction", "rotation_contraction", "projection", "fold", "collapse", "rigid", "identity"]
    stats = {"instances": 0, "failures": 0, "cells": 0, "max_cells": 0, "max_bits": 0}
    bad = []
    with tempfile.TemporaryDirectory() as tmp:
        for it in range(N):
            dom = paper(rng)
            n = rng.randint(1, 6)
            pts = sources(rng, dom, n)
            kind = rng.choice(kinds)
            tgs = targets(rng, pts, kind)
            order = list(range(n))
            rng.shuffle(order)
            c, s = rng.choice(PYTH)
            g = [[c, -s], [s, c]] if rng.random() < 0.5 else [[c, s], [s, -c]]
            L = admissible_L(rng, dom, pts, pts[order[0]])
            st = lambda p: [str(p[0]), str(p[1])]  # noqa: E731
            inst = {"domain": [st(p) for p in dom], "sources": [st(p) for p in pts], "targets": [st(p) for p in tgs],
                    "order": order, "g": [[str(x) for x in r] for r in g], "L": str(L)}
            ip, cp = Path(tmp) / "in.json", Path(tmp) / "cert.json"
            ip.write_text(json.dumps(inst))
            r = subprocess.run([exe, "construct", str(ip), str(cp)], capture_output=True, text=True)
            ok = r.returncode == 0
            if ok:
                summ = json.loads(r.stdout)
                ok = subprocess.run([exe, "check", str(cp)], capture_output=True).returncode == 0
                cert, spec = ref_check.load(cp), ref_check.load(ip)
                ok = ok and not ref_check.check(cert, spec)[1]
                try:
                    ok = ok and brehm_eval.faithful(cert, instance=spec)["agree"]
                except RuntimeError as e:  # the evaluator could not find an exit point
                    ok, r.stderr = False, f"evaluator: {e}"
                stats["cells"] += summ["display_cells"]
                stats["max_cells"] = max(stats["max_cells"], summ["display_cells"])
                stats["max_bits"] = max(stats["max_bits"], summ["max_bits"])
            stats["instances"] += 1
            if not ok:
                stats["failures"] += 1
                bad.append({"instance": inst, "kind": kind, "stderr": r.stderr})
                print("FAIL", json.dumps(bad[-1]), flush=True)
    print(json.dumps({"seed": seed, **stats}))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.setrecursionlimit(100000)
    sys.set_int_max_str_digits(0)
    sys.exit(main())
