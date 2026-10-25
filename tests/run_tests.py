"""Regression and validation suite for brehm_cpp.

    python tests/run_tests.py [path/to/brehm_cpp]

For every instance in examples/: the constructor's exit code must match "expect" (ok -> 0,
invalid -> 3). Every certificate it writes must pass

  * brehm_cpp check (boundary cancellation),
  * tools/ref_check.py --instance (slab decomposition + pairwise continuity, R1-R6 of docs/SPEC.md),
  * tools/brehm_eval.py (agreement with Brehm's map, evaluated without cells, at every cell vertex,
    for the pairs and choices the instance asks for; skipped above EVAL_MAX_CELLS cells),
  * and, where the repository's Python package is present, python/core/exact_certify.py and, for
    instances without construction options, equality with python/core/original_brehm.py.

regression_line and regression_fold must reproduce examples/<name>_cert.json exactly (pieces, domain, data).
Finally every checker must reject a set of corrupted certificates and accept harmless rewrites, and
the evaluator must tell Brehm's map apart from a valid extension that differs from it only by a
hidden fold (SPEC §6: validity alone does not imply faithfulness; vertex agreement does).
"""
import copy
import json
import subprocess
import sys
import tempfile
import time
from fractions import Fraction as F
from pathlib import Path

CPP = Path(__file__).resolve().parents[1]
# python/core sits beside this package in a release, at the repository root in the development tree
PY = next((p / "python" for p in (CPP, *CPP.parents[1:3]) if (p / "python" / "core").is_dir()), None)
sys.path.insert(0, str(CPP / "tools"))
import brehm_eval  # noqa: E402
import ref_check  # noqa: E402

EVAL_MAX_CELLS = 200
try:
    if PY is None:
        raise ImportError("python/core not found")
    sys.path.insert(0, str(PY))
    from core.exact_certify import verify_certificate  # noqa: E402
    from core.original_brehm import certificate_from_json, certificate_to_json, construct as py_construct  # noqa: E402
except ImportError:
    verify_certificate = None

failures = []


def expect(cond, what):
    if not cond:
        failures.append(what)
        print("FAIL", what, flush=True)


def checkers(cert_path, instance_path=None):
    """Verdicts of the three (or two) certificate checkers."""
    out = {"cpp": subprocess.run([EXE, "check", str(cert_path)], capture_output=True).returncode == 0}
    rep, fails = ref_check.check(ref_check.load(cert_path), ref_check.load(instance_path) if instance_path else None)
    out["ref"] = not fails
    if verify_certificate:
        try:
            # read through ref_check's reader, so all three checkers see the same content (SPEC §3)
            out["py"] = verify_certificate(certificate_from_json(ref_check.load(cert_path)))["all_pass"]
        except (ValueError, KeyError, ArithmeticError):
            out["py"] = False
    return out


def run_examples(tmp):
    for inst in sorted((CPP / "examples").glob("*.json")):
        if inst.name.endswith("_cert.json"):
            continue
        spec = json.loads(inst.read_text())
        cert = Path(tmp) / f"{inst.stem}_cert.json"
        t0 = time.perf_counter()
        r = subprocess.run([EXE, "construct", str(inst), str(cert)], capture_output=True, text=True)
        want = {"ok": 0, "invalid": 3}[spec["expect"]]
        expect(r.returncode == want, f"{inst.stem}: exit {r.returncode}, expected {want} ({r.stderr.strip()})")
        if r.returncode != 0:
            print(f"{inst.stem:40s} rejected as expected: {r.stderr.strip()}")
            continue
        summary = json.loads(r.stdout)
        verdicts = checkers(cert, inst)
        expect(all(verdicts.values()), f"{inst.stem}: checkers {verdicts}")
        raw = ref_check.load(cert)
        ev = "skipped"
        if summary["display_cells"] <= EVAL_MAX_CELLS:
            ev = brehm_eval.faithful(raw, instance=ref_check.load(inst))["agree"]
            expect(ev, f"{inst.stem}: certificate differs from Brehm's map")
        same = "n/a"
        if verify_certificate and not any(k in spec for k in ("order", "g", "L")):
            ref = certificate_to_json(py_construct(spec["domain"], spec["sources"], spec["targets"]))
            mine = json.loads(cert.read_text())
            same = all(mine[k] == ref[k] for k in ("domain", "sources", "targets", "pieces"))
            expect(same, f"{inst.stem}: differs from the Python reference constructor")
        print(f"{inst.stem:40s} cells {summary['display_cells']:4d}  checkers {verdicts}  evaluator {ev}  "
              f"python-identical {same}  {time.perf_counter() - t0:.1f}s", flush=True)
        if inst.stem in REGRESSIONS:
            regression(inst.stem, cert)


REGRESSIONS = {"regression_line", "regression_fold"}


def regression(name, cert):
    mine = json.loads(cert.read_text())
    stored = json.loads((CPP / "examples" / f"{name}_cert.json").read_text())
    expect(all(mine[k] == stored[k] for k in ("domain", "sources", "targets", "pieces")),
           f"{name}: certificate differs from the stored regression certificate")


def mutations(cert):
    """(name, certificate, should_pass)."""
    P = cert["pieces"]
    half = lambda s: str(F(s) + F(1, 1000))  # noqa: E731
    yield "nudge a vertex", _with(cert, lambda c: c["pieces"][3]["polygon"][0].__setitem__(0, half(c["pieces"][3]["polygon"][0][0]))), False
    yield "shift one translation", _with(cert, lambda c: c["pieces"][5]["b"].__setitem__(1, half(c["pieces"][5]["b"][1]))), False
    yield "delete a cell", _with(cert, lambda c: c["pieces"].pop(7)), False
    yield "duplicate a cell", _with(cert, lambda c: c["pieces"].append(copy.deepcopy(c["pieces"][2]))), False
    yield "swap two maps", _with(cert, lambda c: _swap(c["pieces"][0], c["pieces"][1])), False
    yield "wrong target", _with(cert, lambda c: c["targets"][1].__setitem__(0, half(c["targets"][1][0]))), False
    yield "non-orthogonal matrix", _with(cert, lambda c: c["pieces"][4]["A"][0].__setitem__(0, str(F(c["pieces"][4]["A"][0][0]) * 2))), False
    yield "move the domain", _with(cert, lambda c: c["domain"][0].__setitem__(0, str(F(c["domain"][0][0]) - 1))), False
    yield "move a source off its cell's image", _with(cert, lambda c: c["sources"][2].__setitem__(1, half(c["sources"][2][1]))), False
    yield "reverse every polygon", _with(cert, lambda c: [p["polygon"].reverse() for p in c["pieces"]]), True
    yield "rotate vertex lists and shuffle cells", _with(cert, lambda c: (c["pieces"].reverse(), [p["polygon"].append(p["polygon"].pop(0)) for p in c["pieces"]])), True
    yield "split a cell in two (same map)", _with(cert, _split), True
    if len(P) > 1:
        yield "split a cell, perturb one half's map", _with(cert, lambda c: (_split(c), c["pieces"][-1]["b"].__setitem__(0, half(c["pieces"][-1]["b"][0])))), False


def _with(cert, f):
    c = copy.deepcopy(cert)
    f(c)
    return c


def _swap(p, q):
    p["A"], q["A"] = q["A"], p["A"]
    p["b"], q["b"] = q["b"], p["b"]


def _split(c):
    """Cut cell 0 along the segment from its vertex 0 to its vertex 2 (convex, at least 3 vertices)."""
    p = c["pieces"][0]
    poly = p["polygon"]
    first, second = poly[0:3], [poly[0]] + poly[2:]
    if len(second) < 3:  # triangle: cut from vertex 0 to the midpoint of the opposite edge
        mid = [str((F(poly[1][i]) + F(poly[2][i])) / 2) for i in range(2)]
        first, second = [poly[0], poly[1], mid], [poly[0], mid, poly[2]]
    p["polygon"] = first
    c["pieces"].append({"name": "split", "polygon": second, "A": copy.deepcopy(p["A"]), "b": copy.deepcopy(p["b"])})


def _fold_corner(c):
    """Fold the tip of a cell at a corner of the domain over a short chord. The result is still a
    valid extension (continuous, piecewise congruent, sources untouched) but not Brehm's map.
    Returns False if no cell has a domain corner whose two cell edges lie on the domain boundary."""
    pt = lambda v: (F(v[0]), F(v[1]))  # noqa: E731
    dom = [pt(v) for v in c["domain"]]
    srcs = [pt(v) for v in c["sources"]]
    on_line = lambda p, q, r: (q[0] - p[0]) * (r[1] - p[1]) == (q[1] - p[1]) * (r[0] - p[0])  # noqa: E731
    for piece in c["pieces"]:
        poly = [pt(v) for v in piece["polygon"]]
        for i, v in enumerate(poly):
            if v not in dom:
                continue
            u, w = poly[i - 1], poly[(i + 1) % len(poly)]
            j = dom.index(v)
            du, dw = dom[j - 1], dom[(j + 1) % len(dom)]
            if not ((on_line(v, du, u) and on_line(v, dw, w)) or (on_line(v, du, w) and on_line(v, dw, u))):
                continue
            eps = F(1, 64)
            p1 = (v[0] + eps * (u[0] - v[0]), v[1] + eps * (u[1] - v[1]))
            p2 = (v[0] + eps * (w[0] - v[0]), v[1] + eps * (w[1] - v[1]))
            tip = [p1, v, p2]
            if any(ref_check.inside(ref_check.convex_ccw(tip), a) for a in srcs):
                continue
            d = (p2[0] - p1[0], p2[1] - p1[1])
            nn = d[0] ** 2 + d[1] ** 2
            R = ((d[0] ** 2 - d[1] ** 2) / nn, 2 * d[0] * d[1] / nn), (2 * d[0] * d[1] / nn, (d[1] ** 2 - d[0] ** 2) / nn)
            A = [[F(x) for x in row] for row in piece["A"]]
            b = [F(x) for x in piece["b"]]
            AR = [[A[r][0] * R[0][k] + A[r][1] * R[1][k] for k in range(2)] for r in range(2)]
            off = (p1[0] - R[0][0] * p1[0] - R[0][1] * p1[1], p1[1] - R[1][0] * p1[0] - R[1][1] * p1[1])
            b2 = [A[r][0] * off[0] + A[r][1] * off[1] + b[r] for r in range(2)]
            piece["polygon"] = [[str(x) for x in q] for q in poly[:i] + [p1, p2] + poly[i + 1:]]
            c["pieces"].append({"name": "folded tip", "polygon": [[str(x) for x in q] for q in tip],
                                "A": [[str(x) for x in row] for row in AR], "b": [str(x) for x in b2]})
            return True
    return False


def run_evaluator_attacks(tmp):
    """Validity does not imply faithfulness; vertex agreement does (SPEC E2). A hidden fold at a
    corner leaves a valid certificate that every checker accepts and the evaluator must reject; a
    refinement (one cell cut in two, same map) must be accepted by both."""
    for name in ("regression_line", "ordinary_random8", "degenerate_collapse", "regression_fold"):
        base = json.loads((Path(tmp) / f"{name}_cert.json").read_text())
        folded = copy.deepcopy(base)
        if not _fold_corner(folded):
            expect(False, f"{name}: no corner cell to fold")
            continue
        for what, cert, valid, faithful in (("hidden fold at a corner", folded, True, False),
                                            ("refinement", _with(base, _split), True, True)):
            path = Path(tmp) / "attack.json"
            path.write_text(json.dumps(cert))
            v = checkers(path)
            ev = brehm_eval.faithful(ref_check.load(path))["agree"]
            expect(all(x == valid for x in v.values()) and ev == faithful,
                   f"{name} / {what}: checkers {v}, evaluator {ev}; expected valid {valid}, faithful {faithful}")
        print(f"evaluator attacks on {name}: done", flush=True)


def run_mutations(tmp):
    for name in ("regression_line", "ordinary_random8", "degenerate_collapse"):
        base = json.loads((Path(tmp) / f"{name}_cert.json").read_text())
        for what, cert, ok in mutations(base):
            path = Path(tmp) / "mut.json"
            path.write_text(json.dumps(cert))
            v = checkers(path)
            expect(all(x == ok for x in v.values()), f"{name} / {what}: verdicts {v}, expected {'pass' if ok else 'reject'}")
        print(f"mutations of {name}: done", flush=True)


def run_starting_isometry_tests(tmp):
    """An inadmissible g can agree at vertices while disagreeing in cell interiors."""
    triangle = [["0", "0"], ["3", "0"], ["0", "3"]]
    base = {"format": "brehm-rational-v1", "name": "starting isometry", "domain": triangle, "sources": [["0", "0"]],
            "targets": [["0", "0"]], "construction": {"order": [0], "L": "1"},
            "pieces": [{"name": "cell", "polygon": triangle,
                        "A": [["-1", "0"], ["0", "-1"]], "b": ["0", "0"]}]}
    cp, ip = Path(tmp) / "g_cert.json", Path(tmp) / "g_instance.json"
    bad_matrices = [
        [["3", "0"], ["0", "3"]],
        [["1", "1/1000000000000000000000000000000"], ["0", "1"]],
        [["1", "1"], ["0", "0"]],  # unit columns that are not perpendicular
        [["1", "0"], ["0", "1/2"]],
        [["1", "0"], ["0", "1"], ["0", "0"]],
        [["1", "0", "0"], ["0", "1", "0"]],
        [["1"], ["0", "1"]], [], "10",
    ]
    for G in bad_matrices:
        cert = copy.deepcopy(base)
        cert["construction"]["g"] = G
        instance = {k: cert[k] for k in ("domain", "sources", "targets")}
        instance.update(cert["construction"])
        cp.write_text(json.dumps(cert))
        ip.write_text(json.dumps(instance))
        expect(all(checkers(cp).values()), "g metadata must not change geometric validity")
        for opts in ([], ["--instance", str(ip)], ["--L", "inf"]):
            p = subprocess.run([sys.executable, "-I", str(CPP / "tools/brehm_eval.py"), str(cp), *opts],
                               capture_output=True, text=True)
            expect(p.returncode == 3 and '"agree": true' not in p.stdout,
                   f"invalid g {G!r} {opts}: evaluator exit {p.returncode}")
        # An invalid instance is rejected even when the certificate has valid choices.
        cert["construction"]["g"] = [["1", "0"], ["0", "1"]]
        cp.write_text(json.dumps(cert))
        p = subprocess.run([sys.executable, "-I", str(CPP / "tools/brehm_eval.py"), str(cp),
                            "--instance", str(ip)], capture_output=True)
        expect(p.returncode == 3, f"invalid instance g {G!r}: evaluator exit {p.returncode}")
        p = subprocess.run([EXE, "construct", str(ip), str(Path(tmp) / "bad_g_output.json")], capture_output=True)
        expect(p.returncode == 3, f"invalid g {G!r}: constructor exit {p.returncode}")
    # Both components of O(2), including non-axis-aligned rational matrices and translation.
    for G in ([["1", "0"], ["0", "1"]], [["0", "-1"], ["1", "0"]],
              [["3/5", "-4/5"], ["4/5", "3/5"]],
              [["-1", "0"], ["0", "1"]], [["3/5", "4/5"], ["4/5", "-3/5"]]):
        instance = {"domain": triangle, "sources": [["1/2", "1/3"]],
                    "targets": [["2/7", "4/9"]], "g": G, "L": "4"}
        ip.write_text(json.dumps(instance))
        p = subprocess.run([EXE, "construct", str(ip), str(cp)], capture_output=True)
        expect(p.returncode == 0, f"valid g {G!r}: constructor exit {p.returncode}")
        if p.returncode == 0:
            expect(all(checkers(cp, ip).values()), f"valid g {G!r}: validity")
            p = subprocess.run([sys.executable, "-I", str(CPP / "tools/brehm_eval.py"), str(cp),
                                "--instance", str(ip), "--points", "4"], capture_output=True, text=True)
            expect(p.returncode == 0 and json.loads(p.stdout)["agree"], f"valid g {G!r}: agreement")
    try:
        brehm_eval.Brehm([(F(0), F(0))], [(F(0), F(0))], ((F(3), F(0)), (F(0), F(3))), F(1))
    except ValueError:
        pass
    else:
        expect(False, "direct Brehm API accepted 3I")
    print("starting isometry tests: done", flush=True)


def main():
    global EXE
    EXE = sys.argv[1] if len(sys.argv) > 1 else str(next(p for p in (CPP / "build" / "brehm_cpp.exe",
                                                                       CPP / "build" / "brehm_cpp") if p.exists()))
    EXE = str(Path(EXE).resolve()) if Path(EXE).is_file() else EXE  # Windows cannot start "dir/exe" relative paths
    with tempfile.TemporaryDirectory() as tmp:
        run_examples(tmp)
        run_mutations(tmp)
        run_evaluator_attacks(tmp)
        run_starting_isometry_tests(tmp)
    print("ALL PASSED" if not failures else f"{len(failures)} FAILURES")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.setrecursionlimit(100000)
    sys.set_int_max_str_digits(0)
    sys.exit(main())
