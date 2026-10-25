// brehm_cpp construct INSTANCE.json CERT.json [--steps]   build and save a certificate
// brehm_cpp check CERT.json                               check a certificate independently
// brehm_cpp orders-exhaustive | orders-anneal ...         insertion-order experiments (orders.cpp)
//
// Instances use the format of python/visualize.py --input (domain, sources, targets), with the
// optional construction choices "order" (0-based indices), "g" (2x2 orthogonal matrix) and "L".
// Certificates use brehm-rational-v1, the format of python/core/original_brehm.py, plus an optional
// "construction" record that the checkers ignore. Exit codes: 0 ok, 1 certificate rejected,
// 2 usage, 3 invalid input, 4 internal error, 5 out of memory.
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

#include "brehm.hpp"
#include "check.hpp"
#include "json.hpp"

namespace {

// GMP aborts when an allocation fails; report it as out of memory (exit 5) instead.
[[noreturn]] void out_of_memory() {
    std::fputs("out of memory\n", stderr);
    std::_Exit(5);
}
void* gmp_alloc(size_t n) {
    void* p = std::malloc(n);
    if (!p) out_of_memory();
    return p;
}
void* gmp_realloc(void* p, size_t, size_t n) {
    void* q = std::realloc(p, n);
    if (!q) out_of_memory();
    return q;
}
void gmp_free(void* p, size_t) { std::free(p); }

std::string slurp(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("cannot read " + path);
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
}

Pt point(const Json& v) {
    if (v.kind != Json::Arr || v.arr.size() != 2) throw std::runtime_error("each point must have two coordinates");
    return {to_rational(v.arr[0]), to_rational(v.arr[1])};
}

std::vector<Pt> points(const Json& v) {
    if (v.kind != Json::Arr) throw std::runtime_error("expected a list of points");
    std::vector<Pt> out;
    for (const Json& p : v.arr) out.push_back(point(p));
    return out;
}

Certificate certificate_from(const Json& j) {
    if (!j.has("format") || j.at("format").text != "brehm-rational-v1")
        throw std::runtime_error("unsupported certificate format");
    Certificate c;
    c.name = j.has("name") ? j.at("name").text : "";
    c.domain = points(j.at("domain"));
    c.sources = points(j.at("sources"));
    c.targets = points(j.at("targets"));
    if (j.at("pieces").kind != Json::Arr) throw std::runtime_error("pieces must be a list");
    for (const Json& p : j.at("pieces").arr) {
        std::vector<Pt> rows = points(p.at("A"));
        if (rows.size() != 2) throw std::runtime_error("a planar isometry must have a 2-by-2 matrix");
        std::string name = p.has("name") ? p.at("name").text : "cell " + std::to_string(c.pieces.size());
        c.pieces.push_back({name, points(p.at("polygon")), Mat{{{rows[0].x, rows[0].y}, {rows[1].x, rows[1].y}}},
                            point(p.at("b"))});
    }
    return c;
}

std::string quoted(const std::string& s) {
    std::string out = "\"";
    for (char ch : s) {
        if (ch == '"' || ch == '\\') out += '\\';
        out += ch;
    }
    return out + "\"";
}

void write_point(std::ostream& o, const Pt& p) { o << "[" << quoted(p.x.get_str()) << "," << quoted(p.y.get_str()) << "]"; }

void write_points(std::ostream& o, const std::vector<Pt>& ps) {
    o << "[";
    for (size_t i = 0; i < ps.size(); ++i) {
        if (i) o << ",";
        write_point(o, ps[i]);
    }
    o << "]";
}

void write_certificate(std::ostream& o, const Certificate& c, const Construction* how = nullptr) {
    o << "{\"format\":\"brehm-rational-v1\",\"name\":" << quoted(c.name) << ",\"domain\":";
    write_points(o, c.domain);
    o << ",\"sources\":";
    write_points(o, c.sources);
    o << ",\"targets\":";
    write_points(o, c.targets);
    o << ",\"pieces\":[";
    for (size_t i = 0; i < c.pieces.size(); ++i) {
        const Piece& p = c.pieces[i];
        if (i) o << ",\n";
        o << "{\"name\":" << quoted(p.name) << ",\"polygon\":";
        write_points(o, p.polygon);
        o << ",\"A\":";
        write_points(o, {{p.A[0][0], p.A[0][1]}, {p.A[1][0], p.A[1][1]}});
        o << ",\"b\":";
        write_point(o, p.b);
        o << "}";
    }
    o << "]";
    if (how) {
        o << ",\n\"construction\":{\"order\":[";
        for (size_t i = 0; i < how->order.size(); ++i) o << (i ? "," : "") << how->order[i];
        o << "],\"g\":";
        write_points(o, {{how->g[0][0], how->g[0][1]}, {how->g[1][0], how->g[1][1]}});
        o << ",\"L\":" << quoted(how->half_width.get_str()) << "}";
    }
    o << "}\n";
}

template <class T>
void write_list(std::ostream& o, const std::vector<T>& v) {
    o << "[";
    for (size_t i = 0; i < v.size(); ++i) o << (i ? "," : "") << v[i];
    o << "]";
}

size_t bits(const Q& v) {
    return std::max(mpz_sizeinbase(v.get_num_mpz_t(), 2), mpz_sizeinbase(v.get_den_mpz_t(), 2));
}

size_t max_bits(const Certificate& c) {
    size_t m = 0;
    for (const Piece& p : c.pieces) {
        for (const Pt& v : p.polygon) m = std::max({m, bits(v.x), bits(v.y)});
        for (const auto& row : p.A) m = std::max({m, bits(row[0]), bits(row[1])});
        m = std::max({m, bits(p.b.x), bits(p.b.y)});
    }
    return m;
}

Options options_from(const Json& in, bool steps) {
    Options opt;
    opt.per_step_display = steps;
    if (in.has("order")) {
        if (in.at("order").kind != Json::Arr) throw std::invalid_argument("order must be a list of indices");
        for (const Json& v : in.at("order").arr) {
            Q k = to_rational(v);
            if (k < 0 || k.get_den() != 1 || !k.get_num().fits_ulong_p())
                throw std::invalid_argument("order must contain nonnegative integers");
            opt.order.push_back(k.get_num().get_ui());
        }
    }
    if (in.has("g")) {
        std::vector<Pt> rows = points(in.at("g"));
        if (rows.size() != 2) throw std::invalid_argument("g must be a 2-by-2 matrix");
        opt.g = Mat{{{rows[0].x, rows[0].y}, {rows[1].x, rows[1].y}}};
    }
    if (in.has("L")) {
        opt.half_width = to_rational(in.at("L"));
        if (opt.half_width <= 0) throw std::invalid_argument("L must be positive");
    }
    return opt;
}

double seconds_since(std::chrono::steady_clock::time_point t0) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

}  // namespace

int orders_main(const std::string& mode, const Json& in, int argc, char** argv);

int main(int argc, char** argv) {
    mp_set_memory_functions(gmp_alloc, gmp_realloc, gmp_free);
    try {
        std::string mode = argc > 1 ? argv[1] : "";
        if (mode == "construct" && (argc == 4 || (argc == 5 && std::string(argv[4]) == "--steps"))) {
            Poly domain;
            std::vector<Pt> src, tgt;
            Options opt;
            try {
                Json in = JsonParser(slurp(argv[2])).parse();
                domain = points(in.at("domain"));
                src = points(in.at("sources"));
                tgt = points(in.at("targets"));
                opt = options_from(in, argc == 5);
            } catch (const std::invalid_argument&) {
                throw;
            } catch (const std::exception& e) {
                throw std::invalid_argument(e.what());
            }
            auto t0 = std::chrono::steady_clock::now();
            Construction r = construct(domain, src, tgt, opt);
            double t = seconds_since(t0);
            std::ofstream out(argv[3], std::ios::binary);
            write_certificate(out, r.certificate, &r);
            if (!out) throw std::runtime_error(std::string("cannot write ") + argv[3]);
            std::cout << "{\"display_cells\":" << r.certificate.pieces.size()
                      << ",\"max_bits\":" << max_bits(r.certificate)
                      << ",\"L\":" << quoted(r.half_width.get_str()) << ",\"full_counts\":";
            write_list(std::cout, r.full_counts);
            std::cout << ",\"display_counts\":";
            write_list(std::cout, r.display_counts);
            std::cout << ",\"seconds\":" << t << "}\n";
            return 0;
        }
        if (mode == "check" && argc == 3) {
            Certificate c;
            try {
                c = certificate_from(JsonParser(slurp(argv[2])).parse());
            } catch (const std::exception& e) {
                throw std::invalid_argument(e.what());
            }
            auto t0 = std::chrono::steady_clock::now();
            CheckReport r = check_certificate(c);
            double t = seconds_since(t0);
            std::cout << "{\"all_pass\":" << (r.all_pass ? "true" : "false") << ",\"pieces\":" << r.pieces
                      << ",\"interior_intervals\":" << r.interior_intervals
                      << ",\"boundary_intervals\":" << r.boundary_intervals << ",\"distinct_maps\":" << r.distinct_maps
                      << ",\"isometric_regions\":" << r.isometric_regions << ",\"max_bits\":" << r.max_bits
                      << ",\"seconds\":" << t << ",\"failures\":[";
            for (size_t i = 0; i < r.failures.size(); ++i) std::cout << (i ? "," : "") << quoted(r.failures[i]);
            std::cout << "]}\n";
            return r.all_pass ? 0 : 1;
        }
        if (mode.rfind("orders-", 0) == 0 && argc >= 3)
            return orders_main(mode, JsonParser(slurp(argv[2])).parse(), argc, argv);
        std::cerr << "usage: brehm_cpp construct INSTANCE.json CERT.json [--steps]\n"
                     "       brehm_cpp check CERT.json\n"
                     "       brehm_cpp orders-exhaustive INSTANCE.json\n"
                     "       brehm_cpp orders-anneal INSTANCE.json ITERS SEED\n";
        return 2;
    } catch (const std::invalid_argument& e) {
        std::cerr << "invalid input: " << e.what() << "\n";
        return 3;
    } catch (const std::bad_alloc&) {
        std::cerr << "out of memory\n";
        return 5;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 4;
    }
}
