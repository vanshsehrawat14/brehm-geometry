// Independent exact certificate checker. It shares only the data types with the constructor.
//
// Method (different from python/core/exact_certify.py, which clips every pair of cells):
// every directed cell edge is grouped by its exact supporting line, the reversed domain boundary
// is added, and on every line the edges are split at all endpoints. If every elementary interval
// has net multiplicity zero, the 1-chain sum(boundary of cells) - boundary(domain) vanishes, so for
// almost every x the number of cells containing x equals 1 inside the domain and 0 outside (each
// cell is a positively oriented convex polygon, winding number 1 inside). That is exact coverage
// with disjoint interiors. Each interior interval then has exactly one cell on each side, and the
// two affine maps are compared at both interval endpoints; vertex-only contacts follow from the
// cyclic order of cells around the vertex. Run time O(E log E) in the number E of edges.
#pragma once

#include <string>
#include <vector>

#include "brehm.hpp"

struct CheckReport {
    bool all_pass = false;
    std::vector<std::string> failures;
    size_t pieces = 0;
    size_t interior_intervals = 0;
    size_t boundary_intervals = 0;
    size_t distinct_maps = 0;        // distinct exact affine maps among the cells
    size_t isometric_regions = 0;    // components of cells glued across edges with equal maps
    size_t max_bits = 0;             // largest numerator/denominator bit length in the cells
};

CheckReport check_certificate(const Certificate& cert);
