"""Exact certificate checker for small Brehm/Petrunin-Yashinski artifacts.

This module is intentionally independent of the production construction and of
Shapely/GEOS.  It verifies a finite polygonal certificate using exact rational
arithmetic:

* each affine matrix is orthogonal;
* pieces lie inside the domain, have disjoint interiors, and their areas sum to
  the domain area;
* affine maps agree exactly on every shared edge;
* every listed source point is mapped exactly to its listed target.

The checker is not a second implementation of the construction.  It consumes a
completed certificate and proves that the certificate describes a continuous
piecewise isometry.  See ``docs/SPEC.md`` (sections 4 and 5) for the
conditions and their proof.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Iterable

Q = Fraction
PointQ = tuple[Q, Q]
MatrixQ = tuple[tuple[Q, Q], tuple[Q, Q]]


def q(value: int | str | Q) -> Q:
    return value if isinstance(value, Q) else Q(value)


def pt(x: int | str | Q, y: int | str | Q) -> PointQ:
    return (q(x), q(y))


@dataclass(frozen=True)
class ExactPiece:
    name: str
    polygon: tuple[PointQ, ...]
    A: MatrixQ
    b: PointQ


@dataclass(frozen=True)
class ExactCertificate:
    name: str
    domain: tuple[PointQ, ...]
    sources: tuple[PointQ, ...]
    targets: tuple[PointQ, ...]
    pieces: tuple[ExactPiece, ...]
    note: str = ""


def add(a: PointQ, b: PointQ) -> PointQ:
    return (a[0] + b[0], a[1] + b[1])


def sub(a: PointQ, b: PointQ) -> PointQ:
    return (a[0] - b[0], a[1] - b[1])


def cross(a: PointQ, b: PointQ) -> Q:
    return a[0] * b[1] - a[1] * b[0]


def dot(a: PointQ, b: PointQ) -> Q:
    return a[0] * b[0] + a[1] * b[1]


def mat_vec(A: MatrixQ, x: PointQ) -> PointQ:
    return (
        A[0][0] * x[0] + A[0][1] * x[1],
        A[1][0] * x[0] + A[1][1] * x[1],
    )


def apply(piece: ExactPiece, x: PointQ) -> PointQ:
    return add(mat_vec(piece.A, x), piece.b)


def signed_area2(poly: Iterable[PointQ]) -> Q:
    pts = tuple(poly)
    if len(pts) < 3:
        return Q(0)
    out = Q(0)
    for a, b in zip(pts, pts[1:] + pts[:1]):
        out += cross(a, b)
    return out


def area(poly: Iterable[PointQ]) -> Q:
    a2 = signed_area2(poly)
    return abs(a2) / 2


def _clean(poly: Iterable[PointQ]) -> tuple[PointQ, ...]:
    out: list[PointQ] = []
    for p in poly:
        if not out or p != out[-1]:
            out.append(p)
    if len(out) > 1 and out[0] == out[-1]:
        out.pop()
    return tuple(out)


def _line_intersection(p1: PointQ, p2: PointQ, q1: PointQ, q2: PointQ) -> PointQ:
    r = sub(p2, p1)
    s = sub(q2, q1)
    den = cross(r, s)
    if den == 0:
        raise ZeroDivisionError("parallel clipping lines")
    t = cross(sub(q1, p1), s) / den
    return (p1[0] + t * r[0], p1[1] + t * r[1])


def convex_intersection(subject: Iterable[PointQ], clipper: Iterable[PointQ]) -> tuple[PointQ, ...]:
    """Exact Sutherland-Hodgman clipping for convex polygons."""
    output = _clean(subject)
    clip = _clean(clipper)
    if len(output) < 3 or len(clip) < 3:
        return ()
    orient = Q(1) if signed_area2(clip) >= 0 else Q(-1)

    def inside(p: PointQ, a: PointQ, b: PointQ) -> bool:
        return orient * cross(sub(b, a), sub(p, a)) >= 0

    for a, b in zip(clip, clip[1:] + clip[:1]):
        if not output:
            break
        inp = output
        output_list: list[PointQ] = []
        prev = inp[-1]
        prev_inside = inside(prev, a, b)
        for cur in inp:
            cur_inside = inside(cur, a, b)
            if cur_inside:
                if not prev_inside:
                    output_list.append(_line_intersection(prev, cur, a, b))
                output_list.append(cur)
            elif prev_inside:
                output_list.append(_line_intersection(prev, cur, a, b))
            prev = cur
            prev_inside = cur_inside
        output = _clean(output_list)
    return output


def point_in_convex_polygon(p: PointQ, poly: Iterable[PointQ]) -> bool:
    pts = _clean(poly)
    if len(pts) < 3:
        return False
    orient = Q(1) if signed_area2(pts) >= 0 else Q(-1)
    return all(orient * cross(sub(b, a), sub(p, a)) >= 0 for a, b in zip(pts, pts[1:] + pts[:1]))


def is_convex_polygon(poly: Iterable[PointQ]) -> bool:
    pts = _clean(poly)
    if len(pts) < 3 or signed_area2(pts) == 0:
        return False
    if any(pts[i] == pts[j] for i in range(len(pts)) for j in range(i + 1, len(pts))):
        return False
    orient = Q(1) if signed_area2(pts) >= 0 else Q(-1)
    # Local turn signs accept star polygons.  Convexity requires every
    # vertex to lie in every consistently oriented edge half-plane.
    return all(
        orient * cross(sub(b, a), sub(p, a)) >= 0
        for a, b in zip(pts, pts[1:] + pts[:1])
        for p in pts
    )


def _segment_overlap(a: PointQ, b: PointQ, c: PointQ, d: PointQ) -> tuple[PointQ, PointQ] | None:
    ab = sub(b, a)
    if cross(ab, sub(c, a)) != 0 or cross(ab, sub(d, a)) != 0:
        return None
    use_x = abs(ab[0]) >= abs(ab[1])
    axis = 0 if use_x else 1
    lo1, hi1 = sorted((a[axis], b[axis]))
    lo2, hi2 = sorted((c[axis], d[axis]))
    lo, hi = max(lo1, lo2), min(hi1, hi2)
    if lo >= hi:
        return None

    def at(value: Q) -> PointQ:
        if a[axis] == b[axis]:
            return a
        t = (value - a[axis]) / (b[axis] - a[axis])
        return (a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1]))

    return (at(lo), at(hi))


def shared_edge_segments(p: Iterable[PointQ], qpoly: Iterable[PointQ]) -> list[tuple[PointQ, PointQ]]:
    ppts = _clean(p)
    qpts = _clean(qpoly)
    out: list[tuple[PointQ, PointQ]] = []
    for a, b in zip(ppts, ppts[1:] + ppts[:1]):
        for c, d in zip(qpts, qpts[1:] + qpts[:1]):
            overlap = _segment_overlap(a, b, c, d)
            if overlap is not None and overlap not in out and (overlap[1], overlap[0]) not in out:
                out.append(overlap)
    return out


def _orthogonal(A: MatrixQ) -> bool:
    col0 = (A[0][0], A[1][0])
    col1 = (A[0][1], A[1][1])
    return dot(col0, col0) == 1 and dot(col1, col1) == 1 and dot(col0, col1) == 0


def _scalars(cert: ExactCertificate) -> Iterable[object]:
    for p in (*cert.domain, *cert.sources, *cert.targets):
        yield from p
    for piece in cert.pieces:
        for p in (*piece.polygon, *piece.A, piece.b):
            yield from p


def verify_certificate(cert: ExactCertificate) -> dict[str, object]:
    # Annotations are not enforced; a float that slips in through the direct API
    # can make inexact comparisons pass (e.g. an off-diagonal 1e-200 squares to 0).
    bad = next((x for x in _scalars(cert) if isinstance(x, bool) or not isinstance(x, (int, Q))), None)
    if bad is not None:
        return {
            "name": cert.name,
            "all_pass": False,
            "failures": [f"non-exact scalar {bad!r}: certificates need int or Fraction"],
            "domain_area": None,
            "total_piece_area": None,
            "piece_count": len(cert.pieces),
            "continuity_edges": 0,
            "max_overlap_area": None,
        }
    failures: list[str] = []
    if len(cert.sources) != len(cert.targets):
        failures.append("source/target length mismatch")

    domain_area = area(cert.domain)
    piece_areas = [area(piece.polygon) for piece in cert.pieces]

    if domain_area <= 0:
        failures.append("domain has zero area")
    if not is_convex_polygon(cert.domain):
        failures.append("domain is not a convex polygon")
    for piece, piece_area in zip(cert.pieces, piece_areas):
        if piece_area <= 0:
            failures.append(f"{piece.name}: zero-area piece")
        if not is_convex_polygon(piece.polygon):
            failures.append(f"{piece.name}: piece is not a convex polygon")
        if not _orthogonal(piece.A):
            failures.append(f"{piece.name}: A is not exactly orthogonal")
        clipped = convex_intersection(piece.polygon, cert.domain)
        if area(clipped) != piece_area:
            failures.append(f"{piece.name}: piece is not contained in the domain")

    total_piece_area = sum(piece_areas, Q(0))
    if total_piece_area != domain_area:
        failures.append(
            f"piece areas do not sum to domain area: {total_piece_area} != {domain_area}"
        )

    max_overlap_area = Q(0)
    continuity_edges = 0
    for i, p in enumerate(cert.pieces):
        for j in range(i + 1, len(cert.pieces)):
            qpiece = cert.pieces[j]
            inter = convex_intersection(p.polygon, qpiece.polygon)
            inter_area = area(inter)
            max_overlap_area = max(max_overlap_area, inter_area)
            if inter_area != 0:
                failures.append(f"{p.name}/{qpiece.name}: positive-area overlap {inter_area}")
            for u, v in shared_edge_segments(p.polygon, qpiece.polygon):
                continuity_edges += 1
                if apply(p, u) != apply(qpiece, u) or apply(p, v) != apply(qpiece, v):
                    failures.append(f"{p.name}/{qpiece.name}: discontinuity on shared edge")

    for idx, (source, target) in enumerate(zip(cert.sources, cert.targets)):
        covering = [piece for piece in cert.pieces if point_in_convex_polygon(source, piece.polygon)]
        if not covering:
            failures.append(f"source {idx}: not covered by any piece")
            continue
        for piece in covering:
            if apply(piece, source) != target:
                failures.append(f"source {idx}: {piece.name} maps to {apply(piece, source)}, not {target}")

    return {
        "name": cert.name,
        "all_pass": not failures,
        "failures": failures,
        "domain_area": domain_area,
        "total_piece_area": total_piece_area,
        "piece_count": len(cert.pieces),
        "continuity_edges": continuity_edges,
        "max_overlap_area": max_overlap_area,
    }


def square_identity_certificate() -> ExactCertificate:
    I = ((Q(1), Q(0)), (Q(0), Q(1)))
    square = (pt(-2, -2), pt(2, -2), pt(2, 2), pt(-2, 2))
    return ExactCertificate(
        name="square_identity", domain=square, sources=square, targets=square,
        pieces=(ExactPiece("square", square, I, pt(0, 0)),),
        note="Identity on a square.",
    )


def two_strip_fold_certificate() -> ExactCertificate:
    """The elementary map (x, y) -> (abs(x), y), with a continuous shared edge."""
    I = ((Q(1), Q(0)), (Q(0), Q(1)))
    R = ((Q(-1), Q(0)), (Q(0), Q(1)))
    return ExactCertificate(
        name="two_strip_fold", domain=(pt(-2, -1), pt(3, -1), pt(3, 1), pt(-2, 1)),
        sources=(pt(-2, 0), pt(0, 0), pt(3, 0)), targets=(pt(2, 0), pt(0, 0), pt(3, 0)),
        pieces=(ExactPiece("left", (pt(-2, -1), pt(0, -1), pt(0, 1), pt(-2, 1)), R, pt(0, 0)),
                ExactPiece("right", (pt(0, -1), pt(3, -1), pt(3, 1), pt(0, 1)), I, pt(0, 0))),
        note="Two rectangles meeting along x = 0.",
    )


def all_certificates() -> tuple[ExactCertificate, ...]:
    return (square_identity_certificate(), two_strip_fold_certificate())


def main() -> int:
    ok = True
    for cert in all_certificates():
        report = verify_certificate(cert)
        ok &= bool(report["all_pass"])
        status = "PASS" if report["all_pass"] else "FAIL"
        print(
            f"{cert.name}: {status} pieces={report['piece_count']} "
            f"area={report['total_piece_area']} continuity_edges={report['continuity_edges']}"
        )
        for failure in report["failures"]:
            print(f"  - {failure}")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
