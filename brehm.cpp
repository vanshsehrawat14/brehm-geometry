// Line-by-line port of python/core/original_brehm.py: same predicates, same branch order, same
// output order, so the two certificates can be compared for exact equality. Only the insertion
// loop differs: intermediate states are restricted to the domain only when asked for.
#include "brehm.hpp"

#include <algorithm>
#include <set>
#include <stdexcept>

namespace {

Q dot(const Pt& a, const Pt& b) { return a.x * b.x + a.y * b.y; }
Q cross(const Pt& a, const Pt& b) { return a.x * b.y - a.y * b.x; }
Pt sub(const Pt& a, const Pt& b) { return {a.x - b.x, a.y - b.y}; }
Q abs_q(const Q& v) { return v < 0 ? Q(-v) : v; }

Q area2(const Poly& p) {
    Q s = 0;
    for (size_t i = 0; i < p.size(); ++i) s += cross(p[i], p[(i + 1) % p.size()]);
    return s;
}

Poly tidy(const Poly& poly) {
    Poly out;
    for (const Pt& p : poly)
        if (out.empty() || p != out.back()) out.push_back(p);
    if (out.size() > 1 && out.front() == out.back()) out.pop_back();
    if (out.size() < 3) return {};
    Q a = area2(out);
    if (a == 0) return {};
    if (a < 0) std::reverse(out.begin(), out.end());
    // Collinear intermediate vertices carry no combinatorial information.
    Poly res;
    size_t n = out.size();
    for (size_t i = 0; i < n; ++i)
        if (cross(sub(out[i], out[(i + n - 1) % n]), sub(out[(i + 1) % n], out[i])) != 0) res.push_back(out[i]);
    return res;
}

bool contains(const Poly& poly, const Pt& p) {
    for (size_t i = 0; i < poly.size(); ++i) {
        const Pt& u = poly[i];
        const Pt& v = poly[(i + 1) % poly.size()];
        if (cross(sub(v, u), sub(p, u)) < 0) return false;
    }
    return true;
}

// Convex polygon clipped to normal.x <= offset, without tolerances.
Poly clip(const Poly& poly, const Pt& normal, const Q& offset) {
    if (poly.empty()) return {};
    Poly out;
    Pt prev = poly.back();
    Q vp = dot(normal, prev) - offset;
    for (const Pt& cur : poly) {
        Q vc = dot(normal, cur) - offset;
        if ((vp < 0 && 0 < vc) || (vc < 0 && 0 < vp)) {
            Q r = vp / (vp - vc);
            out.push_back({prev.x + r * (cur.x - prev.x), prev.y + r * (cur.y - prev.y)});
        }
        if (vc <= 0) out.push_back(cur);
        prev = cur;
        vp = vc;
    }
    return tidy(out);
}

// Extreme points (lexicographic) of the polygon's intersection with normal.x = offset.
std::vector<Pt> section(const Poly& poly, const Pt& normal, const Q& offset) {
    std::set<Pt> points;
    for (size_t i = 0; i < poly.size(); ++i) {
        const Pt& p = poly[i];
        const Pt& q = poly[(i + 1) % poly.size()];
        Q vp = dot(normal, p) - offset, vq = dot(normal, q) - offset;
        if (vp == 0) points.insert(p);
        if ((vp < 0 && 0 < vq) || (vq < 0 && 0 < vp)) {
            Q r = vp / (vp - vq);
            points.insert({p.x + r * (q.x - p.x), p.y + r * (q.y - p.y)});
        }
    }
    if (points.size() < 2) return {};
    return {*points.begin(), *points.rbegin()};
}

}  // namespace

namespace detail {

std::vector<Piece> restrict_to(const std::vector<Piece>& pieces, const Poly& domain) {
    std::vector<Piece> out;
    for (const Piece& piece : pieces) {
        Poly poly = piece.polygon;
        for (size_t i = 0; i < domain.size() && !poly.empty(); ++i) {
            const Pt& u = domain[i];
            Pt edge = sub(domain[(i + 1) % domain.size()], u);
            Pt normal{edge.y, -edge.x};
            poly = clip(poly, normal, dot(normal, u));
        }
        if (!poly.empty()) out.push_back({"cell-" + std::to_string(out.size()), poly, piece.A, piece.b});
    }
    return out;
}

Q default_half_width(const Poly& domain, const Pt& a) {
    Q radius = 1;
    for (const Pt& x : domain) radius = std::max(radius, Q(abs_q(x.x - a.x) + abs_q(x.y - a.y)));
    return 2 * radius;
}

std::vector<Piece> initial(const Poly& domain, const Pt& a, const Pt& b) {
    return initial(a, b, Mat{{{Q(1), Q(0)}, {Q(0), Q(1)}}}, default_half_width(domain, a));
}

// f_1(x) = b + g(phi(x - a)) with phi folding each coordinate at +-L, on the square a + [-3L, 3L]^2.
std::vector<Piece> initial(const Pt& a, const Pt& b, const Mat& g, const Q& L) {
    // lower endpoint, upper endpoint, slope, intercept in centered coordinates
    const std::array<std::array<Q, 4>, 3> intervals{{{Q(-3 * L), Q(-L), Q(-1), Q(-2 * L)},
                                                     {Q(-L), L, Q(1), Q(0)},
                                                     {L, Q(3 * L), Q(-1), Q(2 * L)}}};
    std::vector<Piece> out;
    for (const auto& ix : intervals)
        for (const auto& iy : intervals) {
            Poly poly{{ix[0], iy[0]}, {ix[1], iy[0]}, {ix[1], iy[1]}, {ix[0], iy[1]}};
            for (Pt& p : poly) p = {p.x + a.x, p.y + a.y};
            Mat A{{{g[0][0] * ix[2], g[0][1] * iy[2]}, {g[1][0] * ix[2], g[1][1] * iy[2]}}};
            Pt d{ix[3] - ix[2] * a.x, iy[3] - iy[2] * a.y};
            Pt t{b.x + g[0][0] * d.x + g[0][1] * d.y, b.y + g[1][0] * d.x + g[1][1] * d.y};
            out.push_back({std::to_string(out.size()), poly, A, t});
        }
    return out;
}

std::vector<Piece> insert(const std::vector<Piece>& pieces, const Pt& a, const Pt& b) {
    std::vector<Piece> out;
    out.reserve(2 * pieces.size());
    for (const Piece& piece : pieces) {
        const Mat& Qm = piece.A;
        const Pt& t = piece.b;
        Pt delta = sub(b, t);
        Pt q{Qm[0][0] * delta.x + Qm[1][0] * delta.y, Qm[0][1] * delta.x + Qm[1][1] * delta.y};
        Pt u = sub(q, a);
        Q norm2 = dot(u, u);
        if (norm2 == 0) { out.push_back(piece); continue; }
        Q c = (dot(q, q) - dot(a, a)) / 2;
        bool good = true;
        for (const Pt& p : piece.polygon)
            if (dot(u, p) < c) { good = false; break; }
        if (good) { out.push_back(piece); continue; }
        Poly kept = clip(piece.polygon, {-u.x, -u.y}, -c);
        if (!kept.empty()) out.push_back({"kept", kept, Qm, t});
        std::vector<Pt> edge = section(piece.polygon, u, c);
        if (!edge.empty()) {
            Poly triangle = tidy({a, edge[0], edge[1]});
            if (triangle.empty()) throw std::logic_error("nonzero bisector produced a degenerate cone");
            Mat S{{{Q(1) - 2 * u.x * u.x / norm2, Q(-2 * u.x * u.y / norm2)},
                   {Q(-2 * u.y * u.x / norm2), Q(1) - 2 * u.y * u.y / norm2}}};
            Pt shift{2 * c * u.x / norm2, 2 * c * u.y / norm2};
            Mat newQ;
            for (int i = 0; i < 2; ++i)
                for (int j = 0; j < 2; ++j) newQ[i][j] = Qm[i][0] * S[0][j] + Qm[i][1] * S[1][j];
            Pt newt{t.x + Qm[0][0] * shift.x + Qm[0][1] * shift.y, t.y + Qm[1][0] * shift.x + Qm[1][1] * shift.y};
            out.push_back({"cone", triangle, newQ, newt});
        }
    }
    for (size_t i = 0; i < out.size(); ++i) out[i].name = "cell-" + std::to_string(i);
    return out;
}

void validate(Poly& domain, const std::vector<Pt>& sources, const std::vector<Pt>& targets) {
    if (domain.size() > 1 && domain.front() == domain.back()) domain.pop_back();
    std::set<Pt> distinct(domain.begin(), domain.end());
    if (domain.size() < 3 || distinct.size() != domain.size() || area2(domain) == 0)
        throw std::invalid_argument("domain must be a nondegenerate convex polygon");
    if (area2(domain) < 0) std::reverse(domain.begin(), domain.end());
    for (const Pt& p : domain)
        if (!contains(domain, p)) throw std::invalid_argument("domain must be a simple convex polygon in boundary order");
    domain = tidy(domain);
    if (sources.empty() || sources.size() != targets.size())
        throw std::invalid_argument("need equally many sources and targets, at least one");
    for (const Pt& p : sources)
        if (!contains(domain, p)) throw std::invalid_argument("all sources must be in the domain");
    for (size_t i = 0; i < sources.size(); ++i)
        for (size_t j = 0; j < i; ++j) {
            Pt da = sub(sources[i], sources[j]), db = sub(targets[i], targets[j]);
            if (dot(db, db) > dot(da, da))
                throw std::invalid_argument("expanding prescribed pair " + std::to_string(j) + ", " + std::to_string(i));
        }
}

}  // namespace detail

using namespace detail;

Construction construct(Poly domain, std::vector<Pt> sources, std::vector<Pt> targets, const Options& options) {
    validate(domain, sources, targets);
    const size_t n = sources.size();
    Construction result;
    result.order = options.order;
    if (result.order.empty())
        for (size_t i = 0; i < n; ++i) result.order.push_back(i);
    std::vector<bool> seen(n);
    if (result.order.size() != n) throw std::invalid_argument("order must list every source exactly once");
    for (size_t i : result.order) {
        if (i >= n || seen[i]) throw std::invalid_argument("order must list every source exactly once");
        seen[i] = true;
    }
    const Mat& g = options.g;
    if (g[0][0] * g[0][0] + g[1][0] * g[1][0] != 1 || g[0][1] * g[0][1] + g[1][1] * g[1][1] != 1 ||
        g[0][0] * g[0][1] + g[1][0] * g[1][1] != 0)
        throw std::invalid_argument("g must be exactly orthogonal");
    const Pt& a1 = sources[result.order[0]];
    Q L = options.half_width == 0 ? default_half_width(domain, a1) : options.half_width;
    if (L <= 0) throw std::invalid_argument("L must be positive");
    // Brehm's (II) needs d(a_i, a_1) <= r = L / sqrt(2); the cells must also cover the domain.
    for (const Pt& p : sources)
        if (2 * dot(sub(p, a1), sub(p, a1)) > L * L)
            throw std::invalid_argument("L too small: need L^2 >= 2 |a_i - a_1|^2 for every source");
    for (const Pt& p : domain)
        if (abs_q(p.x - a1.x) > 3 * L || abs_q(p.y - a1.y) > 3 * L)
            throw std::invalid_argument("L too small: the domain must lie in a_1 + [-3L, 3L]^2");
    result.g = g;
    result.half_width = L;
    std::vector<Piece> pieces;
    for (size_t k = 0; k < n; ++k) {
        size_t i = result.order[k];
        pieces = k ? insert(pieces, sources[i], targets[i]) : initial(a1, targets[i], g, L);
        result.full_counts.push_back(pieces.size());
        if (options.per_step_display) result.display_counts.push_back(restrict_to(pieces, domain).size());
    }
    Certificate& cert = result.certificate;
    cert.name = "brehm-step-" + std::to_string(sources.size());
    cert.domain = domain;
    cert.sources = sources;
    cert.targets = targets;
    cert.pieces = restrict_to(pieces, domain);
    return result;
}
