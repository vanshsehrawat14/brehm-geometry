#include "check.hpp"

#include <algorithm>
#include <numeric>
#include <set>
#include <tuple>

namespace {

Q turn(const Pt& o, const Pt& a, const Pt& b) { return (a.x - o.x) * (b.y - o.y) - (a.y - o.y) * (b.x - o.x); }

Q twice_area(const Poly& p) {
    Q s = 0;
    for (size_t i = 0; i < p.size(); ++i) {
        const Pt& u = p[i];
        const Pt& v = p[(i + 1) % p.size()];
        s += u.x * v.y - v.x * u.y;
    }
    return s;
}

// Orients the polygon counter-clockwise; false unless it is a weakly convex polygon with
// positive area and no repeated vertex.
bool orient_convex(Poly& p) {
    if (p.size() < 3) return false;
    if (std::set<Pt>(p.begin(), p.end()).size() != p.size()) return false;
    Q a = twice_area(p);
    if (a == 0) return false;
    if (a < 0) std::reverse(p.begin(), p.end());
    for (size_t i = 0; i < p.size(); ++i)
        for (const Pt& w : p)
            if (turn(p[i], p[(i + 1) % p.size()], w) < 0) return false;
    return true;
}

bool inside(const Poly& ccw, const Pt& x) {
    for (size_t i = 0; i < ccw.size(); ++i)
        if (turn(ccw[i], ccw[(i + 1) % ccw.size()], x) < 0) return false;
    return true;
}

Pt apply(const Piece& c, const Pt& x) {
    return {c.A[0][0] * x.x + c.A[0][1] * x.y + c.b.x, c.A[1][0] * x.x + c.A[1][1] * x.y + c.b.y};
}

bool same_map(const Piece& u, const Piece& v) {
    return u.A[0][0] == v.A[0][0] && u.A[0][1] == v.A[0][1] && u.A[1][0] == v.A[1][0] &&
           u.A[1][1] == v.A[1][1] && u.b == v.b;
}

// Directed edge on its canonical supporting line. kind 1: x + k1*y = k2, parameter y;
// kind 0 (horizontal): y = k2, parameter x.
struct Rec {
    int kind;
    Q k1, k2, lo, hi;
    int sign;
    long owner;  // cell index, or -1 for the reversed domain boundary
};

Rec record(const Pt& p, const Pt& q, long owner, int orientation) {
    Rec r;
    r.owner = owner;
    Q dx = q.x - p.x, dy = q.y - p.y;
    if (dy != 0) {
        r.kind = 1;
        r.k1 = -dx / dy;
        r.k2 = p.x + r.k1 * p.y;
        r.lo = std::min(p.y, q.y);
        r.hi = std::max(p.y, q.y);
        r.sign = (dy > 0 ? 1 : -1) * orientation;
    } else {
        r.kind = 0;
        r.k1 = 0;
        r.k2 = p.y;
        r.lo = std::min(p.x, q.x);
        r.hi = std::max(p.x, q.x);
        r.sign = (dx > 0 ? 1 : -1) * orientation;
    }
    return r;
}

Pt at(const Rec& line, const Q& s) {
    if (line.kind == 1) return {line.k2 - line.k1 * s, s};
    return {s, line.k2};
}

size_t bits(const Q& v) {
    return std::max(mpz_sizeinbase(v.get_num_mpz_t(), 2), mpz_sizeinbase(v.get_den_mpz_t(), 2));
}

struct Failures {
    std::vector<std::string>& list;
    size_t count = 0;
    void add(const std::string& s) {
        if (++count <= 25) list.push_back(s);
    }
};

size_t find(std::vector<size_t>& parent, size_t i) {
    while (parent[i] != i) i = parent[i] = parent[parent[i]];
    return i;
}

}  // namespace

CheckReport check_certificate(const Certificate& cert) {
    CheckReport rep;
    Failures fail{rep.failures};
    rep.pieces = cert.pieces.size();
    if (cert.sources.size() != cert.targets.size()) fail.add("source/target length mismatch");

    Poly domain = cert.domain;
    if (!orient_convex(domain)) {
        fail.add("domain is not a convex polygon with positive area");
        rep.all_pass = false;
        return rep;
    }

    std::vector<Poly> polys;
    polys.reserve(cert.pieces.size());
    bool shapes_ok = true;
    for (size_t i = 0; i < cert.pieces.size(); ++i) {
        const Piece& c = cert.pieces[i];
        Poly p = c.polygon;
        if (!orient_convex(p)) {
            fail.add(c.name + ": not a convex polygon with positive area");
            shapes_ok = false;
        }
        polys.push_back(p);
        const Mat& A = c.A;
        if (A[0][0] * A[0][0] + A[1][0] * A[1][0] != 1 || A[0][1] * A[0][1] + A[1][1] * A[1][1] != 1 ||
            A[0][0] * A[0][1] + A[1][0] * A[1][1] != 0)
            fail.add(c.name + ": linear part is not exactly orthogonal");
        for (const Pt& v : c.polygon) rep.max_bits = std::max({rep.max_bits, bits(v.x), bits(v.y)});
        for (const auto& row : A) rep.max_bits = std::max({rep.max_bits, bits(row[0]), bits(row[1])});
        rep.max_bits = std::max({rep.max_bits, bits(c.b.x), bits(c.b.y)});
    }
    if (!shapes_ok) {
        rep.all_pass = false;
        return rep;
    }

    // Boundary cancellation and continuity.
    std::vector<Rec> recs;
    for (size_t i = 0; i < polys.size(); ++i)
        for (size_t j = 0; j < polys[i].size(); ++j)
            recs.push_back(record(polys[i][j], polys[i][(j + 1) % polys[i].size()], static_cast<long>(i), 1));
    for (size_t j = 0; j < domain.size(); ++j) recs.push_back(record(domain[j], domain[(j + 1) % domain.size()], -1, -1));
    std::sort(recs.begin(), recs.end(), [](const Rec& a, const Rec& b) {
        return std::tie(a.kind, a.k1, a.k2) < std::tie(b.kind, b.k1, b.k2);
    });

    std::vector<size_t> parent(cert.pieces.size());
    std::iota(parent.begin(), parent.end(), 0);
    for (size_t g = 0; g < recs.size();) {
        size_t h = g;
        while (h < recs.size() && recs[h].kind == recs[g].kind && recs[h].k1 == recs[g].k1 && recs[h].k2 == recs[g].k2) ++h;
        std::vector<size_t> by_lo(h - g), by_hi(h - g);
        std::iota(by_lo.begin(), by_lo.end(), g);
        std::iota(by_hi.begin(), by_hi.end(), g);
        std::sort(by_lo.begin(), by_lo.end(), [&](size_t a, size_t b) { return recs[a].lo < recs[b].lo; });
        std::sort(by_hi.begin(), by_hi.end(), [&](size_t a, size_t b) { return recs[a].hi < recs[b].hi; });
        std::vector<Q> cuts;
        for (size_t k = g; k < h; ++k) {
            cuts.push_back(recs[k].lo);
            cuts.push_back(recs[k].hi);
        }
        std::sort(cuts.begin(), cuts.end());
        cuts.erase(std::unique(cuts.begin(), cuts.end()), cuts.end());
        std::vector<size_t> active;
        size_t ia = 0, ib = 0;
        for (size_t k = 0; k + 1 < cuts.size(); ++k) {
            while (ia < by_lo.size() && recs[by_lo[ia]].lo <= cuts[k]) active.push_back(by_lo[ia++]);
            while (ib < by_hi.size() && recs[by_hi[ib]].hi <= cuts[k]) {
                active.erase(std::find(active.begin(), active.end(), by_hi[ib]));
                ++ib;
            }
            if (active.empty()) continue;
            int net = 0;
            for (size_t r : active) net += recs[r].sign;
            Pt p0 = at(recs[g], cuts[k]), p1 = at(recs[g], cuts[k + 1]);
            if (net != 0 || active.size() != 2) {
                fail.add("coverage: multiplicity " + std::to_string(net) + " with " + std::to_string(active.size()) +
                         " edges near (" + p0.x.get_str() + ", " + p0.y.get_str() + ")");
                continue;
            }
            long u = recs[active[0]].owner, v = recs[active[1]].owner;
            if (u < 0 || v < 0) {
                ++rep.boundary_intervals;
                continue;
            }
            ++rep.interior_intervals;
            const Piece& cu = cert.pieces[u];
            const Piece& cv = cert.pieces[v];
            if (apply(cu, p0) != apply(cv, p0) || apply(cu, p1) != apply(cv, p1))
                fail.add(cu.name + "/" + cv.name + ": maps disagree on a shared edge");
            else if (same_map(cu, cv))
                parent[find(parent, u)] = find(parent, v);
        }
        g = h;
    }

    // Prescribed images: every closed cell containing a source must send it to its target.
    std::vector<std::array<Q, 4>> box;
    box.reserve(polys.size());
    for (const Poly& p : polys) {
        std::array<Q, 4> bb{p[0].x, p[0].x, p[0].y, p[0].y};
        for (const Pt& v : p) {
            bb[0] = std::min(bb[0], v.x);
            bb[1] = std::max(bb[1], v.x);
            bb[2] = std::min(bb[2], v.y);
            bb[3] = std::max(bb[3], v.y);
        }
        box.push_back(bb);
    }
    for (size_t s = 0; s < cert.sources.size() && s < cert.targets.size(); ++s) {
        const Pt& a = cert.sources[s];
        if (!inside(domain, a)) fail.add("source " + std::to_string(s) + ": outside the domain");
        size_t covering = 0;
        for (size_t i = 0; i < polys.size(); ++i) {
            if (a.x < box[i][0] || a.x > box[i][1] || a.y < box[i][2] || a.y > box[i][3]) continue;
            if (!inside(polys[i], a)) continue;
            ++covering;
            if (apply(cert.pieces[i], a) != cert.targets[s])
                fail.add("source " + std::to_string(s) + ": " + cert.pieces[i].name + " misses its target");
        }
        if (!covering) fail.add("source " + std::to_string(s) + ": not covered");
    }

    std::set<std::array<Q, 6>> maps;
    for (const Piece& c : cert.pieces) maps.insert({c.A[0][0], c.A[0][1], c.A[1][0], c.A[1][1], c.b.x, c.b.y});
    rep.distinct_maps = maps.size();
    for (size_t i = 0; i < parent.size(); ++i)
        if (find(parent, i) == i) ++rep.isometric_regions;
    if (fail.count > rep.failures.size())
        rep.failures.push_back("... " + std::to_string(fail.count - rep.failures.size()) + " more failures");
    rep.all_pass = fail.count == 0;
    return rep;
}
