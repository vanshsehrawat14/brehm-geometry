// Insertion-order experiments: the display cell count of Brehm's construction depends on the order
// in which the prescribed pairs are inserted (and the first pair also fixes the protected square).
//
//   brehm_cpp orders-exhaustive INSTANCE.json         exact minimum over all n! orders
//   brehm_cpp orders-anneal INSTANCE.json ITERS SEED  simulated annealing over orders
//
// The exhaustive search prunes a prefix as soon as its display count reaches the best complete
// count. That is valid only when counts never decrease along an insertion sequence, which holds for
// collapse data with sources in the open domain: no cell is ever deleted (every crossed cell keeps
// a positive-area part) and every new cone has its apex at the new source. The mode refuses other
// data rather than prune unsoundly.
#include <cmath>
#include <iostream>
#include <limits>
#include <random>

#include "brehm.hpp"
#include "json.hpp"

namespace {

Pt jpoint(const Json& v) { return {to_rational(v.arr.at(0)), to_rational(v.arr.at(1))}; }

struct Instance {
    Poly domain;
    std::vector<Pt> src, tgt;
};

Instance load(const Json& in) {
    Instance I;
    for (const Json& p : in.at("domain").arr) I.domain.push_back(jpoint(p));
    for (const Json& p : in.at("sources").arr) I.src.push_back(jpoint(p));
    for (const Json& p : in.at("targets").arr) I.tgt.push_back(jpoint(p));
    detail::validate(I.domain, I.src, I.tgt);
    return I;
}

size_t display(const std::vector<Piece>& pieces, const Poly& domain) {
    return detail::restrict_to(pieces, domain).size();
}

size_t evaluate(const Instance& I, const std::vector<size_t>& order) {
    std::vector<Piece> pieces = detail::initial(I.domain, I.src[order[0]], I.tgt[order[0]]);
    for (size_t k = 1; k < order.size(); ++k) pieces = detail::insert(pieces, I.src[order[k]], I.tgt[order[k]]);
    return display(pieces, I.domain);
}

bool strictly_inside(const Poly& ccw, const Pt& x) {
    for (size_t i = 0; i < ccw.size(); ++i) {
        const Pt& u = ccw[i];
        const Pt& v = ccw[(i + 1) % ccw.size()];
        if ((v.x - u.x) * (x.y - u.y) - (v.y - u.y) * (x.x - u.x) <= 0) return false;
    }
    return true;
}

struct Exhaustive {
    const Instance& I;
    size_t best = std::numeric_limits<size_t>::max();
    std::vector<size_t> best_order{}, order{};
    std::vector<bool> used{};
    size_t nodes = 0;

    void dfs(const std::vector<Piece>& pieces) {
        size_t n = I.src.size();
        if (order.size() == n) {
            size_t c = display(pieces, I.domain);
            if (c < best) {
                best = c;
                best_order = order;
            }
            return;
        }
        for (size_t i = 0; i < n; ++i) {
            if (used[i]) continue;
            std::vector<Piece> next = order.empty() ? detail::initial(I.domain, I.src[i], I.tgt[i])
                                                    : detail::insert(pieces, I.src[i], I.tgt[i]);
            ++nodes;
            if (display(next, I.domain) >= best) continue;
            used[i] = true;
            order.push_back(i);
            dfs(next);
            order.pop_back();
            used[i] = false;
        }
    }
};

void print_order(const std::vector<size_t>& o) {
    std::cout << "[";
    for (size_t i = 0; i < o.size(); ++i) std::cout << (i ? "," : "") << o[i];
    std::cout << "]";
}

}  // namespace

int orders_main(const std::string& mode, const Json& in, int argc, char** argv) {
    Instance I = load(in);
    size_t n = I.src.size();
    std::vector<size_t> natural(n), reversed(n);
    for (size_t i = 0; i < n; ++i) natural[i] = i, reversed[i] = n - 1 - i;
    size_t c_nat = evaluate(I, natural), c_rev = evaluate(I, reversed);
    if (mode == "orders-exhaustive") {
        for (size_t i = 1; i < n; ++i)
            if (I.tgt[i] != I.tgt[0]) throw std::invalid_argument("exhaustive pruning needs collapse data");
        for (const Pt& a : I.src)
            if (!strictly_inside(I.domain, a)) throw std::invalid_argument("exhaustive pruning needs interior sources");
        Exhaustive E{I};
        E.used.assign(n, false);
        // Start from the better of the two fixed orders; the search only looks for strictly better.
        E.best = std::min(c_nat, c_rev);
        E.best_order = c_nat <= c_rev ? natural : reversed;
        E.dfs({});
        std::cout << "{\"natural\":" << c_nat << ",\"reversed\":" << c_rev << ",\"minimum\":" << E.best
                  << ",\"nodes\":" << E.nodes << ",\"argmin\":";
        print_order(E.best_order);
        std::cout << "}\n";
        return 0;
    }
    if (mode == "orders-anneal" && argc == 5) {
        long iters = std::stol(argv[3]);
        std::mt19937_64 rng(std::stoull(argv[4]));
        std::vector<size_t> cur = c_nat <= c_rev ? natural : reversed, best = cur;
        double cur_c = static_cast<double>(std::min(c_nat, c_rev)), best_c = cur_c;
        double t0 = 0.05 * cur_c, t1 = 0.002 * cur_c;
        std::uniform_int_distribution<size_t> pick(0, n - 1);
        std::uniform_real_distribution<double> unit(0.0, 1.0);
        for (long it = 0; it < iters; ++it) {
            double T = t0 * std::pow(t1 / t0, static_cast<double>(it) / static_cast<double>(iters));
            std::vector<size_t> cand = cur;
            size_t i = pick(rng), j = pick(rng);
            if (i == j) continue;
            if (unit(rng) < 0.5) {
                std::swap(cand[i], cand[j]);
            } else {
                size_t v = cand[i];
                cand.erase(cand.begin() + static_cast<long>(i));
                cand.insert(cand.begin() + static_cast<long>(j), v);
            }
            double c = static_cast<double>(evaluate(I, cand));
            if (c <= cur_c || unit(rng) < std::exp((cur_c - c) / T)) {
                cur = cand;
                cur_c = c;
                if (c < best_c) {
                    best_c = c;
                    best = cand;
                }
            }
        }
        std::cout << "{\"natural\":" << c_nat << ",\"reversed\":" << c_rev << ",\"best\":" << best_c << ",\"order\":";
        print_order(best);
        std::cout << "}\n";
        return 0;
    }
    std::cerr << "usage: brehm_cpp orders-exhaustive INSTANCE.json | orders-anneal INSTANCE.json ITERS SEED\n";
    return 2;
}
