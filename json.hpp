// Minimal JSON reader for the instance and certificate formats, with exact numbers.
// Numbers keep their literal text so "0.1" becomes exactly 1/10, as in the Python reader.
// Duplicate keys and nesting deeper than 64 levels are rejected, so every reader of a file sees
// the same content (Python's json keeps the last duplicate, a naive reader the first).
#pragma once

#include <gmpxx.h>

#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

struct Json {
    enum Kind { Null, Bool, Num, Str, Arr, Obj } kind = Null;
    std::string text;  // number literal or string contents
    bool flag = false;
    std::vector<Json> arr;
    std::vector<std::pair<std::string, Json>> obj;

    const Json& at(const std::string& key) const {
        for (const auto& kv : obj)
            if (kv.first == key) return kv.second;
        throw std::runtime_error("missing key: " + key);
    }
    bool has(const std::string& key) const {
        for (const auto& kv : obj)
            if (kv.first == key) return true;
        return false;
    }
};

inline bool is_digit(char c) { return c >= '0' && c <= '9'; }

class JsonParser {
public:
    explicit JsonParser(const std::string& s) : s_(s) {}
    Json parse() {
        Json v = value();
        ws();
        if (i_ != s_.size()) fail("trailing characters");
        return v;
    }

private:
    const std::string& s_;
    size_t i_ = 0;
    int depth_ = 0;

    struct Nest {
        int& depth;
        explicit Nest(int& d) : depth(d) { ++depth; }
        ~Nest() { --depth; }
    };

    [[noreturn]] void fail(const std::string& what) {
        throw std::runtime_error("JSON: " + what + " at offset " + std::to_string(i_));
    }
    void ws() {
        while (i_ < s_.size() && (s_[i_] == ' ' || s_[i_] == '\n' || s_[i_] == '\r' || s_[i_] == '\t')) ++i_;
    }
    char peek() {
        ws();
        if (i_ >= s_.size()) fail("unexpected end");
        return s_[i_];
    }
    void expect(char c) {
        if (peek() != c) fail(std::string("expected '") + c + "'");
        ++i_;
    }
    bool literal(const char* word) {
        size_t n = std::char_traits<char>::length(word);
        if (s_.compare(i_, n, word) == 0) { i_ += n; return true; }
        return false;
    }
    std::string string() {
        expect('"');
        std::string out;
        while (true) {
            if (i_ >= s_.size()) fail("unterminated string");
            char c = s_[i_++];
            if (c == '"') return out;
            if (c != '\\') { out += c; continue; }
            if (i_ >= s_.size()) fail("bad escape");
            char e = s_[i_++];
            switch (e) {
                case '"': case '\\': case '/': out += e; break;
                case 'b': out += '\b'; break;
                case 'f': out += '\f'; break;
                case 'n': out += '\n'; break;
                case 'r': out += '\r'; break;
                case 't': out += '\t'; break;
                case 'u': {
                    if (i_ + 4 > s_.size()) fail("bad \\u escape");
                    unsigned code = 0;
                    for (size_t k = 0; k < 4; ++k, ++i_) {
                        char h = s_[i_];
                        int d = is_digit(h) ? h - '0' : (h >= 'a' && h <= 'f') ? h - 'a' + 10 : (h >= 'A' && h <= 'F') ? h - 'A' + 10 : -1;
                        if (d < 0) fail("bad \\u escape");
                        code = 16 * code + static_cast<unsigned>(d);
                    }
                    if (code > 0x7f) fail("non-ASCII \\u escape is not supported");
                    out += static_cast<char>(code);
                    break;
                }
                default: fail("bad escape");
            }
        }
    }
    Json value() {
        Nest nest(depth_);
        if (depth_ > 64) fail("nesting deeper than 64 levels");
        Json v;
        char c = peek();
        if (c == '{') {
            ++i_;
            v.kind = Json::Obj;
            if (peek() == '}') { ++i_; return v; }
            std::set<std::string> keys;
            while (true) {
                std::string key = string();
                if (!keys.insert(key).second) fail("duplicate key \"" + key + "\"");
                expect(':');
                v.obj.emplace_back(key, value());
                if (peek() == ',') { ++i_; continue; }
                expect('}');
                return v;
            }
        }
        if (c == '[') {
            ++i_;
            v.kind = Json::Arr;
            if (peek() == ']') { ++i_; return v; }
            while (true) {
                v.arr.push_back(value());
                if (peek() == ',') { ++i_; continue; }
                expect(']');
                return v;
            }
        }
        if (c == '"') { v.kind = Json::Str; v.text = string(); return v; }
        if (literal("true")) { v.kind = Json::Bool; v.flag = true; return v; }
        if (literal("false")) { v.kind = Json::Bool; return v; }
        if (literal("null")) return v;
        size_t start = i_;
        while (i_ < s_.size() && std::string("+-0123456789.eE").find(s_[i_]) != std::string::npos) ++i_;
        if (start == i_) fail("unexpected character");
        v.kind = Json::Num;
        v.text = s_.substr(start, i_ - start);
        return v;
    }
};

// Exact value of a decimal literal such as "-12", "0.125" or "1e-3". The explicit exponent is
// limited to |e| <= 1000 so that a short literal cannot demand a huge number.
inline mpq_class parse_decimal(const std::string& t) {
    size_t i = 0;
    bool neg = false;
    if (i < t.size() && (t[i] == '+' || t[i] == '-')) neg = (t[i++] == '-');
    std::string digits;
    long exponent = 0;
    bool any = false;
    while (i < t.size() && is_digit(t[i])) { digits += t[i++]; any = true; }
    if (i < t.size() && t[i] == '.') {
        ++i;
        while (i < t.size() && is_digit(t[i])) { digits += t[i++]; --exponent; any = true; }
    }
    if (!any) throw std::runtime_error("not a number: " + t);
    if (i < t.size() && (t[i] == 'e' || t[i] == 'E')) {
        ++i;
        bool eneg = false;
        if (i < t.size() && (t[i] == '+' || t[i] == '-')) eneg = (t[i++] == '-');
        if (i == t.size()) throw std::runtime_error("not a number: " + t);
        long e = 0;
        for (; i < t.size() && is_digit(t[i]); ++i) {
            e = 10 * e + (t[i] - '0');
            if (e > 1000) throw std::runtime_error("exponent out of range (|e| <= 1000): " + t);
        }
        exponent += eneg ? -e : e;
    }
    if (i != t.size()) throw std::runtime_error("not a number: " + t);
    mpz_class num(digits.empty() ? "0" : digits, 10), scale;
    mpz_ui_pow_ui(scale.get_mpz_t(), 10, static_cast<unsigned long>(exponent < 0 ? -exponent : exponent));
    mpq_class q = exponent < 0 ? mpq_class(num, scale) : mpq_class(num * scale);
    q.canonicalize();
    return neg ? mpq_class(-q) : q;
}

// A coordinate: a JSON number, or a string holding a decimal or "p/q" (p an optionally signed
// integer, q a positive integer). No spaces. Booleans are rejected. docs/SPEC.md §3 has the grammar.
inline mpq_class to_rational(const Json& v) {
    if (v.kind == Json::Num) return parse_decimal(v.text);
    if (v.kind != Json::Str) throw std::runtime_error("coordinates must be numbers or rational strings");
    const std::string& t = v.text;
    size_t slash = t.find('/');
    if (slash == std::string::npos) return parse_decimal(t);
    size_t start = (t[0] == '+' || t[0] == '-') ? 1 : 0;
    bool ok = start < slash && slash + 1 < t.size();
    for (size_t i = start; ok && i < t.size(); ++i) ok = i == slash || is_digit(t[i]);
    if (!ok) throw std::runtime_error("not a number: " + t);
    mpz_class num(t.substr(t[0] == '+' ? 1 : 0, slash - (t[0] == '+' ? 1 : 0)), 10), den(t.substr(slash + 1), 10);
    if (den == 0) throw std::runtime_error("zero denominator: " + v.text);
    mpq_class q(num, den);
    q.canonicalize();
    return q;
}
