"""Rational planar reference solver for Brehm's 1981 construction.

The folded outer square is part of the algorithm, not a plotting window.
See docs/SPEC.md, section 2. No geometry from the certificate checker is
used here; only the certificate data types are
shared. This is an implementation of known mathematics, not a new theorem.
"""

from fractions import Fraction as F
from itertools import product

from .exact_certify import ExactCertificate, ExactPiece


def _point(p):
    if len(p) != 2:
        raise ValueError("each point must have two coordinates")
    if any(isinstance(x, (float, bool)) for x in p):
        raise ValueError("use integers, rational strings, or Fraction, not floats")
    return tuple(F(x) for x in p)


def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1]


def _sub(a, b):
    return a[0] - b[0], a[1] - b[1]


def _cross(a, b):
    return a[0] * b[1] - a[1] * b[0]


def _area2(poly):
    return sum((_cross(p, q) for p, q in zip(poly, poly[1:] + poly[:1])), F(0))


def _tidy(poly):
    out = []
    for p in poly:
        if not out or p != out[-1]:
            out.append(p)
    if len(out) > 1 and out[0] == out[-1]:
        out.pop()
    if len(out) < 3 or _area2(out) == 0:
        return ()
    if _area2(out) < 0:
        out.reverse()
    # Collinear intermediate vertices carry no combinatorial information.
    out = [p for i, p in enumerate(out)
           if _cross(_sub(p, out[i-1]), _sub(out[(i+1) % len(out)], p))]
    return tuple(out)


def _contains(poly, p):
    return all(_cross(_sub(v, u), _sub(p, u)) >= 0
               for u, v in zip(poly, poly[1:] + poly[:1]))


def _clip(poly, normal, offset):
    """Convex polygon clipped to normal.x <= offset, without tolerances."""
    out = []
    if not poly:
        return ()
    prev = poly[-1]
    vp = _dot(normal, prev) - offset
    for cur in poly:
        vc = _dot(normal, cur) - offset
        if (vp < 0 < vc) or (vc < 0 < vp):
            ratio = vp / (vp - vc)
            out.append(tuple(prev[j] + ratio * (cur[j] - prev[j]) for j in range(2)))
        if vc <= 0:
            out.append(cur)
        prev, vp = cur, vc
    return _tidy(out)


def _section(poly, normal, offset):
    points = set()
    for p, q in zip(poly, poly[1:] + poly[:1]):
        vp, vq = _dot(normal, p) - offset, _dot(normal, q) - offset
        if vp == 0:
            points.add(p)
        if (vp < 0 < vq) or (vq < 0 < vp):
            ratio = vp / (vp - vq)
            points.add(tuple(p[j] + ratio * (q[j] - p[j]) for j in range(2)))
    points = sorted(points)
    return (points[0], points[-1]) if len(points) >= 2 else ()


def _restrict(pieces, domain):
    out = []
    for piece in pieces:
        poly = piece.polygon
        for u, v in zip(domain, domain[1:] + domain[:1]):
            edge = _sub(v, u)
            normal = edge[1], -edge[0]
            poly = _clip(poly, normal, _dot(normal, u))
            if not poly:
                break
        if poly:
            out.append(ExactPiece(f"cell-{len(out)}", poly, piece.A, piece.b))
    return tuple(out)


def _validate(domain, sources, targets):
    domain = tuple(_point(p) for p in domain)
    if len(domain) > 1 and domain[0] == domain[-1]:
        domain = domain[:-1]
    if len(domain) < 3 or len(set(domain)) != len(domain) or not _area2(domain):
        raise ValueError("domain must be a nondegenerate convex polygon")
    if _area2(domain) < 0:
        domain = tuple(reversed(domain))
    if not all(_contains(domain, p) for p in domain):
        raise ValueError("domain must be a simple convex polygon in boundary order")
    domain = _tidy(domain)
    sources, targets = tuple(map(_point, sources)), tuple(map(_point, targets))
    if not sources or len(sources) != len(targets):
        raise ValueError("need equally many sources and targets, at least one")
    if not all(_contains(domain, p) for p in sources):
        raise ValueError("all sources must be in the domain")
    for i, a in enumerate(sources):
        for j in range(i):
            da, db = _sub(a, sources[j]), _sub(targets[i], targets[j])
            if _dot(db, db) > _dot(da, da):
                raise ValueError(f"expanding prescribed pair {j}, {i}")
    return domain, sources, targets


def _initial(domain, a, b):
    radius = max(F(1), *(sum(abs(x[j] - a[j]) for j in range(2)) for x in domain))
    length = 2 * radius
    # lower endpoint, upper endpoint, slope, intercept in centered coordinates
    intervals = ((-3*length, -length, F(-1), -2*length),
                 (-length, length, F(1), F(0)),
                 (length, 3*length, F(-1), 2*length))
    out = []
    for ix, iy in product(intervals, repeat=2):
        x0, x1, sx, tx = ix
        y0, y1, sy, ty = iy
        poly = tuple((x+a[0], y+a[1]) for x, y in
                     ((x0, y0), (x1, y0), (x1, y1), (x0, y1)))
        out.append(ExactPiece(str(len(out)), poly, ((sx, F(0)), (F(0), sy)),
                              (b[0]+tx-sx*a[0], b[1]+ty-sy*a[1])))
    return tuple(out)


def _insert(pieces, a, b):
    out = []
    for piece in pieces:
        Q, t, poly = piece.A, piece.b, piece.polygon
        delta = _sub(b, t)
        q = tuple(sum(Q[j][i] * delta[j] for j in range(2)) for i in range(2))
        u = _sub(q, a)
        norm2 = _dot(u, u)
        if norm2 == 0:
            out.append(piece)
            continue
        c = (_dot(q, q) - _dot(a, a)) / 2
        if all(_dot(u, p) >= c for p in poly):
            out.append(piece)
            continue
        kept = _clip(poly, (-u[0], -u[1]), -c)
        if kept:
            out.append(ExactPiece("kept", kept, Q, t))
        edge = _section(poly, u, c)
        if edge:
            triangle = _tidy((a, *edge))
            if not triangle:
                raise ArithmeticError("nonzero bisector produced a degenerate cone")
            S = tuple(tuple(F(i == j) - 2*u[i]*u[j]/norm2 for j in range(2))
                      for i in range(2))
            shift = tuple(2*c*u[i]/norm2 for i in range(2))
            newQ = tuple(tuple(sum(Q[i][k]*S[k][j] for k in range(2))
                               for j in range(2)) for i in range(2))
            newt = tuple(t[i] + sum(Q[i][k]*shift[k] for k in range(2)) for i in range(2))
            out.append(ExactPiece("cone", triangle, newQ, newt))
    return tuple(ExactPiece(f"cell-{i}", p.polygon, p.A, p.b) for i, p in enumerate(out))


def construct_steps(domain, sources, targets, *, name="brehm"):
    """Yield (restricted exact certificate, full-square cell count) per anchor.

    The caller may independently verify each certificate. No finite test of
    this implementation is a formal proof of success on all admissible inputs.
    """
    domain, sources, targets = _validate(domain, sources, targets)
    pieces = _initial(domain, sources[0], targets[0])
    for i in range(len(sources)):
        if i:
            pieces = _insert(pieces, sources[i], targets[i])
        cert = ExactCertificate(f"{name}-step-{i+1}", domain, sources[:i+1], targets[:i+1],
                                _restrict(pieces, domain))
        yield cert, len(pieces)


def construct(domain, sources, targets, *, name="brehm"):
    """Return the final rational certificate; verify it with exact_certify."""
    result = None
    for result, _ in construct_steps(domain, sources, targets, name=name):
        pass
    return result


def certificate_to_json(cert):
    def point(p):
        return [str(x) for x in p]
    return {
        "format": "brehm-rational-v1", "name": cert.name,
        "domain": list(map(point, cert.domain)),
        "sources": list(map(point, cert.sources)),
        "targets": list(map(point, cert.targets)),
        "pieces": [{"name": p.name, "polygon": list(map(point, p.polygon)),
                    "A": list(map(point, p.A)), "b": point(p.b)} for p in cert.pieces],
    }


def certificate_from_json(raw):
    """Load an untrusted witness. Call verify_certificate before using it."""
    if raw.get("format") != "brehm-rational-v1":
        raise ValueError("unsupported certificate format")
    def matrix(rows):
        parsed = tuple(map(_point, rows))
        if len(parsed) != 2:
            raise ValueError("a planar isometry must have a 2-by-2 matrix")
        return parsed
    return ExactCertificate(str(raw["name"]), tuple(map(_point, raw["domain"])),
                            tuple(map(_point, raw["sources"])), tuple(map(_point, raw["targets"])),
                            tuple(ExactPiece(str(p["name"]), tuple(map(_point, p["polygon"])),
                                             matrix(p["A"]), _point(p["b"]))
                                  for p in raw["pieces"]))
