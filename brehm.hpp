// Exact planar Brehm extension (C++/GMP port of python/core/original_brehm.py).
#pragma once

#include <gmpxx.h>

#include <array>
#include <string>
#include <vector>

using Q = mpq_class;

struct Pt {
    Q x, y;
};
inline bool operator==(const Pt& a, const Pt& b) { return a.x == b.x && a.y == b.y; }
inline bool operator!=(const Pt& a, const Pt& b) { return !(a == b); }
inline bool operator<(const Pt& a, const Pt& b) { return a.x < b.x || (a.x == b.x && a.y < b.y); }

using Poly = std::vector<Pt>;
using Mat = std::array<std::array<Q, 2>, 2>;

// One convex cell carrying the isometry x -> A x + b.
struct Piece {
    std::string name;
    Poly polygon;
    Mat A;
    Pt b;
};

struct Certificate {
    std::string name;
    Poly domain;
    std::vector<Pt> sources, targets;
    std::vector<Piece> pieces;
};

// The choices Brehm's proof leaves free (Brehm 1981, pp. 188-189). The defaults reproduce
// python/core/original_brehm.py exactly.
struct Options {
    std::vector<size_t> order;  // insertion order as indices into the sources; empty: as given
    Mat g{{{Q(1), Q(0)}, {Q(0), Q(1)}}};  // linear part of the initial isometry, exactly orthogonal
    Q half_width = 0;           // L: folds at +-L about the first inserted source, so Brehm's r = L/sqrt(2);
                                // 0 picks L = 2R, R = max(1, L1 radius of the domain about that source)
    bool per_step_display = false;  // also restrict every intermediate state to the domain
};

struct Construction {
    Certificate certificate;            // final cells restricted to the domain
    std::vector<size_t> full_counts;    // protected-square cells after each insertion
    std::vector<size_t> display_counts; // cells meeting the domain after each insertion (if requested)
    std::vector<size_t> order;          // the insertion order used
    Mat g;                              // the linear part of g used
    Q half_width;                       // the L used
};

// Same contract and the same cell order as original_brehm.construct; throws std::invalid_argument
// on inadmissible input or options.
Construction construct(Poly domain, std::vector<Pt> sources, std::vector<Pt> targets, const Options& options = {});

// The steps of construct, for drivers that try many insertion orders.
namespace detail {
void validate(Poly& domain, const std::vector<Pt>& sources, const std::vector<Pt>& targets);
Q default_half_width(const Poly& domain, const Pt& a);
std::vector<Piece> initial(const Poly& domain, const Pt& a, const Pt& b);
std::vector<Piece> initial(const Pt& a, const Pt& b, const Mat& g, const Q& half_width);
std::vector<Piece> insert(const std::vector<Piece>& pieces, const Pt& a, const Pt& b);
std::vector<Piece> restrict_to(const std::vector<Piece>& pieces, const Poly& domain);
}  // namespace detail
