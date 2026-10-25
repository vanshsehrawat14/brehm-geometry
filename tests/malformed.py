"""Malformed and hostile input: every tool must fail explicitly, quickly and without a crash.

    python tests/malformed.py [path/to/brehm_cpp]

Runs brehm_cpp construct / check, tools/ref_check.py and tools/brehm_eval.py on broken instances
and certificates (bad JSON, duplicate keys, 100,000 nested brackets, numbers outside the grammar,
exponents that would need gigabytes, bad options, wrong shapes) and on one legitimate instance with
5,000-digit coordinates. Each run must end within TIMEOUT seconds with an expected exit code
(3 invalid input, 2 usage, 1 rejected certificate) and must not print a Python traceback.
"""
import json
import subprocess
import sys
import tempfile
import time
from pathlib import Path

CPP = Path(__file__).resolve().parents[1]
TOOLS = CPP / "tools"
TIMEOUT = 60
failures = []

SQUARE = [["-1", "-1"], ["1", "-1"], ["1", "1"], ["-1", "1"]]
GOOD = {"domain": SQUARE, "sources": [["0", "0"], ["1/2", "0"]], "targets": [["0", "0"], ["1/4", "0"]]}


def instance(**changes):
    out = dict(GOOD)
    out.update(changes)
    return json.dumps(out)


INSTANCES = {  # name: file text; construct must exit 3
    "not JSON": "domain: square",
    "empty file": "",
    "top level is a list": "[1, 2]",
    "missing domain": json.dumps({k: v for k, v in GOOD.items() if k != "domain"}),
    "duplicate key": instance()[:-1] + ', "sources": [["1", "1"]]}',
    "100,000 nested brackets": instance(note="x")[:-1] + ', "deep": ' + "[" * 100000 + "]" * 100000 + "}",
    "65 nesting levels (limit 64)": instance(note="x")[:-1] + ', "deep": ' + "[" * 64 + "]" * 64 + "}",
    "exponent 10^999999999 in a string": instance(sources=[["1e999999999", "0"], ["1/2", "0"]]),
    "exponent 10^999999999 as a JSON number": instance()[:-1].replace('["0", "0"], ["1/2"', '[1e999999999, 0], ["1/2"', 1) + "}",
    "zero denominator": instance(sources=[["1/0", "0"], ["1/2", "0"]]),
    "spaces inside a number": instance(sources=[["1 / 2", "0"], ["1/2", "0"]]),
    "negative denominator": instance(sources=[["1/-2", "0"], ["1/2", "0"]]),
    "underscore digits": instance(sources=[["1_0", "0"], ["1/2", "0"]]),
    "NaN literal": instance()[:-1].replace('["0", "0"], ["1/2"', '[NaN, 0], ["1/2"', 1) + "}",
    "boolean coordinate": instance(sources=[[True, "0"], ["1/2", "0"]]),
    "three coordinates": instance(sources=[["0", "0", "0"], ["1/2", "0"]]),
    "null domain": instance(domain=None),
    "fewer targets than sources": instance(targets=[["0", "0"]]),
    "bad unicode escape": instance(note="x").replace('"x"', '"\\uZZZZ"'),
    "truncated file": instance()[: len(instance()) // 2],
    "trailing characters": instance() + " x",
    "order with a fraction": instance(order=["1/2", "0"]),
    "order with a repeat": instance(order=[0, 0]),
    "order too short": instance(order=[0]),
    "g not orthogonal": instance(g=[["1", "1"], ["0", "1"]]),
    "L = 0": instance(L="0"),
    "L negative": instance(L="-3"),
    "L below Brehm's bound": instance(L="1/10"),
    "non-convex paper": instance(domain=[["0", "0"], ["2", "0"], ["1", "1/4"], ["2", "2"], ["0", "2"]]),
    "expanding pair": instance(targets=[["0", "0"], ["2", "0"]]),
}

CERT = {"format": "brehm-rational-v1", "name": "c", "domain": SQUARE, "sources": [["0", "0"]], "targets": [["0", "0"]],
        "pieces": [{"name": "p", "polygon": SQUARE, "A": [["1", "0"], ["0", "1"]], "b": ["0", "0"]}]}


def certificate(**changes):
    out = dict(CERT)
    out.update(changes)
    return json.dumps(out)


CERTIFICATES = {  # name: (file text, exit codes allowed for the validity checkers)
    "valid identity (control)": (certificate(), {0}),
    "not JSON": ("pieces", {3}),
    "duplicate pieces key": (certificate()[:-1] + ', "pieces": []}', {3}),
    "100,000 nested brackets": (certificate()[:-1] + ', "deep": ' + "[" * 100000 + "]" * 100000 + "}", {3}),
    "64 nesting levels (allowed)": (certificate()[:-1] + ', "deep": ' + "[" * 63 + "]" * 63 + "}", {0}),
    "65 nesting levels": (certificate()[:-1] + ', "deep": ' + "[" * 64 + "]" * 64 + "}", {3}),
    "65 levels with a scalar innermost": (certificate()[:-1] + ', "deep": ' + "[" * 63 + "1" + "]" * 63 + "}", {3}),
    "exponent 10^999999999": (certificate(sources=[["1e999999999", "0"]]), {3}),
    "zero denominator": (certificate(targets=[["0", "1/0"]]), {3}),
    "pieces not a list": (certificate(pieces={"p": 1}), {3}),
    "piece without a map": (certificate(pieces=[{"polygon": SQUARE, "A": [["1", "0"], ["0", "1"]]}]), {3}),
    "piece without a name (allowed)": (certificate(pieces=[{"polygon": SQUARE, "A": [["1", "0"], ["0", "1"]], "b": ["0", "0"]}]), {0}),
    "two-point polygon": (certificate(pieces=[dict(CERT["pieces"][0], polygon=SQUARE[:2])]), {1}),
    "matrix with one row": (certificate(pieces=[dict(CERT["pieces"][0], A=[["1", "0"]])]), {1, 3}),
    "wrong format tag": (certificate(format="brehm-float"), {1, 3}),
}


def run(cmd, allowed, what):
    t0 = time.perf_counter()
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=TIMEOUT)
    except subprocess.TimeoutExpired:
        failures.append(f"{what}: no answer within {TIMEOUT} s")
        print("FAIL", failures[-1], flush=True)
        return
    dt = time.perf_counter() - t0
    ok = r.returncode in allowed and "Traceback" not in r.stderr
    if not ok:
        failures.append(f"{what}: exit {r.returncode}, allowed {sorted(allowed)}; {r.stderr.strip()[-300:]}")
    print(f"{'ok  ' if ok else 'FAIL'} {what:70s} exit {r.returncode}  {dt:.2f}s  {r.stderr.strip()[:90]}", flush=True)


def main():
    exe = sys.argv[1] if len(sys.argv) > 1 else str(next(p for p in (CPP / "build" / "brehm_cpp.exe",
                                                                    CPP / "build" / "brehm_cpp") if p.exists()))
    exe = str(Path(exe).resolve()) if Path(exe).is_file() else exe  # Windows cannot start "dir/exe" relative paths
    py = [sys.executable, "-I"]
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        out = tmp / "out.json"
        for name, text in INSTANCES.items():
            path = tmp / "inst.json"
            path.write_text(text, encoding="utf-8")
            run([exe, "construct", str(path), str(out)], {3}, f"construct / {name}")
        path = tmp / "inst.json"
        path.write_text(instance(), encoding="utf-8")
        run([exe, "construct", str(path), str(out), "--fast"], {2}, "construct / unknown flag")
        run([exe, "construct", str(path)], {2}, "construct / missing output path")
        run([exe], {2}, "no command")
        run([exe, "construct", str(tmp / "missing.json"), str(out)], {3}, "construct / missing file")
        for name, (text, allowed) in CERTIFICATES.items():
            path = tmp / "cert.json"
            path.write_text(text, encoding="utf-8")
            run([exe, "check", str(path)], allowed, f"check / {name}")
            run(py + [str(TOOLS / "ref_check.py"), str(path)], allowed, f"ref_check / {name}")
            # the evaluator decides faithfulness, not validity: it may accept a structurally bad file
            run(py + [str(TOOLS / "brehm_eval.py"), str(path)], allowed | {0, 1}, f"brehm_eval / {name}")
        run(py + [str(TOOLS / "ref_check.py")], {2}, "ref_check / no arguments")
        run(py + [str(TOOLS / "brehm_eval.py"), str(tmp / "cert.json"), "--bogus", "1"], {2}, "brehm_eval / unknown option")
        # Legitimate input with 5,000-digit coordinates: everything must accept it (no digit limits).
        big = "1" + "0" * 5000
        path.write_text(json.dumps({"domain": [["-" + big, "-" + big], [big, "-" + big], [big, big], ["-" + big, big]],
                                    "sources": [["0", "0"], [big[:-1], "0"]], "targets": [["0", "0"], ["0", big[:-1]]]}))
        run([exe, "construct", str(path), str(out)], {0}, "construct / 5,000-digit coordinates")
        run([exe, "check", str(out)], {0}, "check / 5,000-digit coordinates")
        run(py + [str(TOOLS / "ref_check.py"), str(out), "--instance", str(path)], {0}, "ref_check / 5,000-digit coordinates")
        run(py + [str(TOOLS / "brehm_eval.py"), str(out), "--instance", str(path)], {0}, "brehm_eval / 5,000-digit coordinates")
    print("ALL PASSED" if not failures else f"{len(failures)} FAILURES")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
