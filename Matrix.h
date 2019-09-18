//
// Created by inFinity on 2019-05-14.
//

#pragma once

#include <vector>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <queue>

#include <eigen3/Eigen/Dense>

#include "bliss/defs.hh"
#include "bliss/graph.hh"
#include "bliss/timer.hh"
#include "bliss/utils.hh"
#include "bliss/bignum.hh"
#include "bliss/uintseqhash.hh"

#include "hash.h"
#include "frac.h"
#include "tree.h"
//#include "BloomFilter.h"

//#define PRINT_DEBUG

using uint = unsigned int;

/// A Symmetric Matrix
class Matrix {
public:
    struct CanonicalGraph {
        uint vtx_num;
        std::vector<int> colors;
        std::vector<unsigned long> adj_bits;

        void init(uint vn) {
            vtx_num = vn;
            colors.resize(vtx_num);
            auto bit_size = vtx_num * (vtx_num - 1) / 2;
            auto ret_size = bit_size / 64;
            ret_size = (bit_size % 64 == 0) ? ret_size : ret_size + 1;
            adj_bits.resize(ret_size);
        }

        auto setAdj(const std::vector<bool> &bits) {
            const int bit_num = 64;
            auto bit_size = bits.size();
            auto ret_size = bit_size / bit_num;
            ret_size = (bit_size % bit_num == 0) ? ret_size : ret_size + 1;
            for (auto i = 0u; i < ret_size; ++i) {
                auto &word = adj_bits[i];
                for (auto j = 0; j < bit_num; ++j) {
                    auto bits_idx = i * bit_num + j;
                    if (bits_idx >= bit_size) break;
                    word = word << 1u;
                    word |= bits[bits_idx] ? 1ul : 0ul;
                }
            }
        }

        CanonicalGraph &operator=(const CanonicalGraph &) = default;

        CanonicalGraph &operator=(CanonicalGraph &&g) noexcept {
            std::swap(vtx_num, g.vtx_num);
            std::swap(colors, g.colors);
            std::swap(adj_bits, g.adj_bits);
            return *this;
        }
    };

public:
    enum HashFunction {
        Eigen, Bliss, BlissHash, TreeHash, EVC, PR, EigenPR, None
    };

    int n;
//        std::vector <std::vector<int>> adj;
    std::vector<int> adj;
    bool is_connected{false};

    // eigen
    std::vector<int> eigen;
    std::vector<int> degree;

    // bliss

    CanonicalGraph cg;
    size_t hash_value;

    // Tree
    Tree t;

    // EVC
    std::vector<double> evc;

    // PageRank
    std::vector<double> pr;
//    std::vector<Fraction> pr;
    size_t pr_hash_value;

    HashFunction hash_func;

    Matrix() = default;

    explicit Matrix(int vertNum, HashFunction hf = None) : n(vertNum), hash_func(hf) {
        auto len = n * (n + 1) / 2;
        adj.resize(len);
        degree.resize(n);
        setHash(hf);
    }

    Matrix(const Matrix &a) {
        n = a.n;
        adj = a.adj;
        hash_func = a.hash_func;
        is_connected = a.is_connected;

        eigen = a.eigen;
        degree = a.degree;

        hash_value = a.hash_value;
        cg = a.cg;

        t = a.t;

        evc = a.evc;

        pr = a.pr;
        pr_hash_value = a.pr_hash_value;
    }

    Matrix(Matrix &a) {
        std::swap(n, a.n);
        std::swap(adj, a.adj);
        std::swap(hash_func, a.hash_func);
        std::swap(is_connected, a.is_connected);

        std::swap(eigen, a.eigen);
        std::swap(degree, a.degree);

        std::swap(t, a.t);

        std::swap(hash_value, a.hash_value);
        cg = a.cg;

        std::swap(evc, a.evc);

        std::swap(pr, a.pr);
        std::swap(pr_hash_value, a.pr_hash_value);
    }

    inline void setHash(HashFunction hf) {
        hash_func = hf;
        hash_value = 0;
        switch (hf) {
            case Eigen:
                eigen.resize(n);
                break;
            case Bliss:
                cg.init(n);
                break;
            case BlissHash:
                break;
            case TreeHash:
                t = Tree(n, adj);
                break;
            case EVC:
                evc.resize(n, 1);
                break;
            case PR:
                pr.resize(n, 1);
                break;
            case EigenPR:
                eigen.resize(n);
                pr.resize(n, 1);
                break;
            default:
                break;
        }
    }

    inline void setDegree() {
        std::fill(degree.begin(), degree.end(), 0);
        for (auto i = 0; i < n; ++i) {
            for (auto j = 0; j < n; ++j) {
                if (get(i, j)) degree[i]++;
            }
        }
    }

    inline void reset() {
        degree.clear();
        eigen.clear();
        evc.clear();
        pr.clear();
    }

    inline int ij2idx(int i, int j) const {
        if (i > j) std::swap(i, j);
        return n * i + j - i * (i + 1) / 2;
    }

    inline int &operator()(int i, int j) {
        return adj[ij2idx(i, j)];
    }

    inline int get(int i, int j) const {
        return adj[ij2idx(i, j)];
    }

    auto track() const {
        int ret = 0;
        for (auto i = 0; i < n; ++i) {
            ret += get(i, i);
        }
        return ret;
    }

    Matrix &operator=(const Matrix &a) {
        n = a.n;
        adj = a.adj;
        hash_func = a.hash_func;
        is_connected = a.is_connected;

        eigen = a.eigen;
        degree = a.degree;

        hash_value = a.hash_value;
        cg = a.cg;

        evc = a.evc;

        pr = a.pr;
        pr_hash_value = a.pr_hash_value;
        return *this;
    }

    Matrix &operator=(Matrix &&a) noexcept {
        std::swap(n, a.n);
        std::swap(adj, a.adj);
        std::swap(hash_func, a.hash_func);
        std::swap(is_connected, a.is_connected);

        std::swap(eigen, a.eigen);
        std::swap(degree, a.degree);

        std::swap(hash_value, a.hash_value);
        cg = a.cg;

        std::swap(evc, a.evc);

        std::swap(pr, a.pr);
        return *this;
    }

    void clone(const Matrix &a) {
        n = a.n;
        std::copy(a.adj.begin(), a.adj.end(), adj.begin());
        hash_func = a.hash_func;
        is_connected = a.is_connected;

        eigen.resize(a.eigen.size());
        std::copy(a.eigen.begin(), a.eigen.end(), eigen.begin());
        degree.resize(a.degree.size());
        std::copy(a.degree.begin(), a.degree.end(), degree.begin());

        hash_value = a.hash_value;

        evc.resize(a.evc.size());
        std::copy(a.evc.begin(), a.evc.end(), evc.begin());

        pr.resize(a.pr.size());
        std::copy(a.pr.begin(), a.pr.end(), pr.begin());
    }


    Matrix operator+(const Matrix &b) {
        assert(n == b.n);
        Matrix ret(n);
        for (auto i = 0; i < n; ++i) {
            for (auto j = i; j < n; ++j) {
                ret(i, j) = get(i, j) + b.get(i, j);
            }
        }
        return std::move(ret);
    }

    Matrix operator*(const Matrix &b) {
        Matrix ret(n);
        for (auto i = 0; i < n; ++i) {
            for (auto j = i; j < n; ++j) {
                auto &target = ret(i, j);
                for (auto it = 0; it < n; ++it) {
                    target += get(i, it) * b.get(it, j);
                }
            }
        }
        return std::move(ret);
    }

    /// Inspired by https://math.stackexchange.com/questions/864604/checking-connectivity-of-adjacency-matrix
    bool isConnected() {
        auto isContainFullRow = [&](Matrix &A) {
            for (auto i = 0; i < A.n; ++i) {
                auto tRet = true;
                for (auto j = 0; j < A.n; ++j) {
                    if (!A.get(i, j)) {
                        tRet = false;
                        break;
                    }
                }
                if (tRet) return true;
            }
            return false;
        };
        Matrix A(n);
        A.clone(*this); // A = M+I
        for (auto i = 0; i < n; ++i) A(i, i) = 1;

        if (isContainFullRow(A)) is_connected = true;
        for (auto i = 0; i < n - 1; ++i) {
            A = multiply(A, *this); // A = (M+I)^(i+2)
            if (isContainFullRow(A)) is_connected = true;
        }
        return is_connected;
    }

    bool isConnectedBFS() const {
        std::vector<bool> visited(n - 1);
        std::queue<uint> frontier;
        frontier.push(0u);
        while (!frontier.empty()) {
            auto x = frontier.front();
            for (auto i = 0u; i < n; ++i) {
                if (get(x, i) && !visited[i]) {
                    frontier.push(i);
                }
            }
            frontier.pop();
        }
        for (auto b: visited) {
            if (!b) return false;
        }
        return true;
    }

    bool isEmpty() const {
        for (auto a: adj) {
            if (a != 0) return true;
        }
        return false;
    }

    inline bool isCanonicalGraphEqual(const CanonicalGraph &rg) const {
        return (cg.vtx_num == rg.vtx_num) && (cg.colors == rg.colors) && (cg.adj_bits == rg.adj_bits);

    }

    template<typename T=int>
    friend Matrix operator*(const Matrix &a, const T &v) {
        auto n = a.n;
        Matrix ret(n);
        for (auto i = 0; i < n; ++i) {
            for (auto j = i; j < n; ++j) {
                ret(i, j) = a.get(i, j) * v;
            }
        }
        return std::move(ret);
    }

    template<typename T=int>
    friend Matrix operator*(const T &v, Matrix &a) {
        return std::move(a * v);
    }

    friend void mPlus(Matrix &a, int v) {
        auto n = a.n;
        for (auto i = 0; i < n; ++i) {
            a(i, i) = a.get(i, i) + v;
        }
    }

    friend Matrix multiply(const Matrix &a, const Matrix &b) {
        auto n = a.n;
        Matrix ret(n);
        for (auto i = 0; i < n; ++i) {
            for (auto j = i; j < n; ++j) {
                auto &target = ret(i, j);
                for (auto it = 0; it < n; ++it) {
                    target += a.get(i, it) * b.get(it, j);
                }
            }
        }
        return std::move(ret);
    }

    friend auto multiply(const Matrix &a, const std::vector<double> &vec) {
        auto n = a.n;
        assert(vec.size() == n);
        std::vector<double> ret(n, 0);
        for (auto i = 0; i < n; ++i) {
            for (auto j = 0; j < n; ++j) {
                ret[i] += a.get(i, j) * vec[j];
            }
        }
        return std::move(ret);
    }

    friend auto multiply(const Matrix &a, const std::vector<Fraction> &vec) {
        auto n = a.n;
        assert(vec.size() == n);
        std::vector<Fraction> ret(n, Fraction(0, 1));
        for (auto i = 0; i < n; ++i) {
            for (auto j = 0; j < n; ++j) {
                ret[i] += a.get(i, j) * vec[j];
            }
        }
        return ret;
    }

    static auto normalize(std::vector<double> &vec) {
        double sum = 0.0f;
        for (auto e: vec) {
            sum += e * e;
        }
        sum = sqrt(sum);
        for (auto &e: vec) {
            e /= sum;
            e *= 10000;
            e = floor(e);
            e /= 10000;
        }
        return sum;
    }

    friend bool operator==(const Matrix &lhs, const Matrix &rhs) {
        bool ret = lhs.hash_func == rhs.hash_func;
        if (!ret) return ret;
        switch (lhs.hash_func) {
            case Matrix::Eigen:
                return lhs.hash_value == rhs.hash_value;
//                if (lhs.hash_value != rhs.hash_value) return false;
//                if (lhs.eigen.empty() || rhs.eigen.empty()) return false;
//                if (lhs.degree.empty() || rhs.degree.empty()) return false;
//                ret = lhs.degree == rhs.degree;
//                if (!ret) return false;
//                return lhs.eigen == rhs.eigen;
            case Matrix::Bliss:
                if (lhs.hash_value != rhs.hash_value) return false;
                return lhs.isCanonicalGraphEqual(rhs.cg);
            case Matrix::BlissHash:
                return lhs.hash_value == rhs.hash_value;
            case Matrix::TreeHash:
                return lhs.hash_value == rhs.hash_value;
            case Matrix::EVC:
                if (lhs.evc.empty() || rhs.evc.empty()) return false;
                if (lhs.degree.empty() || rhs.degree.empty()) return false;
                ret = lhs.degree == rhs.degree;
                if (!ret) return false;
                return lhs.evc == rhs.evc;
            case Matrix::PR:
                if (lhs.pr.empty() || rhs.pr.empty()) return false;
                if (lhs.degree.empty() || rhs.degree.empty()) return false;
                ret = lhs.degree == rhs.degree;
                if (!ret) return false;
                return lhs.pr_hash_value == rhs.pr_hash_value;
            case Matrix::EigenPR:
                if (lhs.degree.empty() || rhs.degree.empty()) return false;
                if (lhs.eigen.empty() || rhs.eigen.empty()) return false;
                ret = lhs.degree == rhs.degree;
                if (!ret) return false;
                ret = lhs.eigen == rhs.eigen;
                if (!ret) return false;
                if (lhs.n >= 8) {
                    return lhs.pr_hash_value == rhs.pr_hash_value;
                } else {
                    return ret;
                }
            default:
                return false;
        }
    }

    friend std::ostream &operator<<(std::ostream &os, const Matrix &m) {
//        os << "[";
        for (auto a: m.adj) {
            os << a;
        }
#ifdef PRINT_DEBUG
        os << "]";
        os << ": [";
        switch (m.hash_func) {
            case Matrix::Eigen:
                for (auto e: m.eigen) {
                    os << " " << e;
                }
                break;
            case Matrix::Bliss:
                os << " " << m.hash_value;
                break;
            case Matrix::TreeHash:
                os << m.t << " " << m.hash_value;
                break;
            case Matrix::EVC:
                for (auto e: m.evc) {
                    os << " " << e;
                }
                break;
            case Matrix::PR:
                for (auto e: m.pr) {
                    os << " " << e;
                }
                break;
            case Matrix::EigenPR:
                for (auto e: m.eigen) {
                    os << " " << e;
                }
                os << " ] [";
                for (auto e: m.pr) {
                    os << " " << e;
                }
                break;
            default:
                break;
        }
        os << " ]";
#endif
        return os;
    }

    void printMatrix() {
        for (auto i = 0; i < n; ++i) {
            for (auto j = 0; j < n; ++j) {
                std::cout << adj[ij2idx(i, j)] << " ";
            }
            std::cout << " --- " << degree[i] << std::endl;
        }
        std::cout << "[";
        switch (hash_func) {
            case Matrix::Eigen:
                for (auto e: eigen) {
                    std::cout << " " << e;
                }
                break;
            case Matrix::Bliss:
                std::cout << " " << hash_value;
                break;
            case Matrix::TreeHash:
                std::cout << t << " " << hash_value;
                break;
            case Matrix::EVC:
                for (auto e: evc) {
                    std::cout << " " << e;
                }
                break;
            case Matrix::PR:
                for (auto e: pr) {
                    std::cout << " " << e;
                }
                break;
            case Matrix::EigenPR:
                for (auto e: eigen) {
                    std::cout << " " << e;
                }
                std::cout << " ] [";
                for (auto e: pr) {
                    std::cout << " " << e;
                }
                break;
            default:
                break;
        }
        std::cout << " ]\n";
    }

    void sortDegree() {
        std::sort(degree.begin(), degree.end());
    }

    /// Characteristic Polynomial
    void charPoly() {
//        using EigenMat = Eigen::MatrixXd;
//        EigenMat mA(n, n);
//        for (auto i = 0u; i < n; ++i) {
//            for (auto j = 0u; j < n; ++j) {
//                mA(i, j) = get(i, j);
//            }
//        }
//        EigenMat mC = mA;
//        for (auto k = 1; k <= n; ++k) {
//            if (k > 1)  {
//                mC = mA * (mC.array() + eigen[n -k +1]).matrix();
//            }
//            eigen[n - k] = -mC.trace()/k;
//        }
        Matrix A(n);
        Matrix C(n);
        A.clone(*this);
        C.clone(*this);
        for (auto k = 1; k <= n; ++k) {
            if (k > 1) {
//                    C = A * (C + c[n - k + 1] * I);
                mPlus(C, eigen[n - k + 1]);
                C.clone(multiply(A, C));
            }
            eigen[n - k] = -C.track() / k;
        }
    }

    void treeHash() {
        // Connected
        // Edge Num == Vertex num -1

        using Ngh = std::vector<uint>;
        std::vector<Ngh> nghs(n);
        for (auto i = 0; i < n - 1; ++i) {
            for (auto j = i + 1; j < n; ++j) {
                if (get(i, j)) {
                    nghs[i].push_back(j);
                    nghs[j].push_back(i);
                }
            }
        }

        auto root = [&]() {
            std::vector<bool> visited(n);
            visited[0] = true;
        };
    }

    void blissing() {
        bliss::AbstractGraph *g = new bliss::Graph(n);
        for (auto i = 0; i < n; ++i) {
            g->change_color(i, 1);
        }

        for (auto i = 0; i < n - 1; ++i) {
            for (auto j = i + 1; j < n; ++j) {
                if (get(i, j)) {
                    g->add_edge(i, j);
                    g->add_edge(j, i);
                }
            }
        }
        bliss::Stats stats;
        auto cl = g->canonical_form(stats, [](void *, const uint, const uint *) {}, stdout);
        auto bg = dynamic_cast<bliss::Graph *>(g->permute(cl));
        hash_value = bg->get_hash();

        auto vn = cg.vtx_num;
        std::vector<bool> bits(vn * (vn - 1) / 2);
        auto setBit = [&](uint i, uint j) {
            if (i > j) std::swap(i, j);
            auto idx = (vn - 1) * i - i * (i + 1) / 2 + j - 1;
            bits[idx] = true;
        };
        for (auto i = 0u; i < cg.vtx_num; ++i) {
            cg.colors[i] = bg->vertices[i].color;
        }
        for (auto v = 0u; v < n; ++v) {
            auto &edges = bg->vertices[v].edges;
            for (auto u: edges) {
                setBit(v, u);
            }
        }
        cg.setAdj(bits);
        delete g;
        delete bg;
    }

    void blissingHash() {
        bliss::AbstractGraph *g = new bliss::Graph(n);
        for (auto i = 0; i < n; ++i) {
            g->change_color(i, 1);
        }

        for (auto i = 0; i < n - 1; ++i) {
            for (auto j = i + 1; j < n; ++j) {
                if (get(i, j)) {
                    g->add_edge(i, j);
                    g->add_edge(j, i);
                }
            }
        }
        bliss::Stats stats;
        auto cl = g->canonical_form(stats, [](void *, const uint, const uint *) {}, stdout);
        auto bg = dynamic_cast<bliss::Graph *>(g->permute(cl));
        hash_value = bg->get_hash();

        delete g;
        delete bg;
    }

    void computeEvc() {
        for (double nv = 0.0;;) {
            evc = multiply(*this, evc);
            auto norm_value = normalize(evc);
            if (abs(nv - norm_value) < 0.0001) break;
            nv = norm_value;
        }
        std::sort(evc.begin(), evc.end());
    }

    void prHash() {
        using std::hash;
#ifdef PRINT_DEBUG
        std::cout << "[";
#endif
        auto bit_cast = [](const auto &p) {
            auto num = *((uint64_t *) (&p));
            const uint64_t mask = 0xfful;
            auto lower = num & mask;
            auto higher = (num & (~mask)) >> 8;
            if (lower > 0x7ful) {
                higher += 1;
            }
#ifdef PRINT_DEBUG
            std::cout << " " << std::hex << higher << std::dec;
#endif
            return higher;
        };
        auto ret = (hash<uint64_t>()(bit_cast(pr.front())) << 7) >> 31;
        for (auto i = 1ul; i < pr.size(); ++i) {
            ret ^= ((hash<uint64_t>()(bit_cast(pr[i]))) << (7 + i)) >> (31 - i);
        }
#ifdef PRINT_DEBUG
        std::cout << " ]\n";
#endif
        pr_hash_value = ret;
    }

    void pagerank() {
        ///
        auto printV = [](const auto &vec) {
            std::cout << "[";
            for (auto &v: vec) {
                std::cout << " " << v;
            }
            std::cout << " ]";
        };
        ///
        auto cmpPr = [&](const auto &pr1, const auto &pr2) {
            double diff = 0.0f;
            for (auto i = 0u; i < pr1.size(); ++i) {
                diff += fabs(pr1[i] - pr2[i]);
            }
//            std::cout << diff << std::endl;
            return (diff < 0.001);
        };
        auto it = 1;
        for (std::vector<double> pr2(n, 0); !cmpPr(pr, pr2); pr.swap(pr2), it++) {
            for (auto i = 0u; i < pr.size(); ++i) {
                pr[i] /= degree[i];
            }
            pr2 = multiply(*this, pr);
            for (auto &p: pr2) {
//                p *= Fraction(17, 20);
//                p += Fraction(3, 20) / n;
                p += 0.85;
                p += 0.15 / n;
            }
#ifdef PRINT_DEBUG
            std::cout << "It[" << it << "]: " << *this << std::endl;
#endif
            if (it >= 5) break;
        }
//        for (auto &p: pr) {
//            p *= 10000;
//            p = round(p);
//            p /= 10000;
//        }
        std::sort(pr.begin(), pr.end());
        prHash();
    }

    void hash() {
        setDegree();
        switch (hash_func) {
            case Matrix::Eigen:
                charPoly();
                sortDegree();
                hash_value = 0;
                for (int i : eigen) {
                    hashUpdate(hash_value, i);
//                    hash_value ^= (std::hash<int>()(eigen[i]) << (17 - i)) >> (11 + i);
                }
                for (int i : degree) {
                    hashUpdate(hash_value, i);
//                    hash_value ^= (std::hash<int>()(degree[i]) << (19 - i)) >> (13 + i);
                }
                break;
            case Matrix::Bliss:
                blissing();
                break;
            case Matrix::BlissHash:
                blissingHash();
                break;
            case Matrix::TreeHash:
                hash_value = t.hash();
                break;
            case Matrix::EVC:
                computeEvc();
                sortDegree();
                break;
            case Matrix::PR:
                pagerank();
                sortDegree();
                break;
            case Matrix::EigenPR:
                charPoly();
                if (n >= 8)
                    pagerank();
                sortDegree();
                break;
            default:
                break;
        }
    }

};

struct MatrixHasher {
    size_t operator()(const Matrix &k) const {
        using std::hash;
        size_t ret;
        switch (k.hash_func) {
            case Matrix::Eigen:
                return k.hash_value;
            case Matrix::Bliss:
                return k.hash_value;
            case Matrix::BlissHash:
                return k.hash_value;
            case Matrix::TreeHash:
                return k.hash_value;
            case Matrix::EVC:
                ret = (hash<double>()(k.evc.front()) << 17) >> 13;
                for (auto i = 1u; i < k.evc.size(); ++i) {
                    ret ^= (hash<double>()(k.evc[i]) << (17 - i)) >> (13 + i);
                }
                for (auto i = 0u; i < k.degree.size(); ++i) {
                    ret ^= (hash<int>()(k.degree[i]) << (19 - i)) >> (13 + i);
                }
                return ret;
            case Matrix::PR:
                return k.pr_hash_value;
            case Matrix::EigenPR:
                ret = (hash<int>()(k.eigen[0]) << 17) >> 11;
                for (int i = 1; i < k.eigen.size(); ++i) {
                    ret ^= (hash<int>()(k.eigen[i]) << (17 - i)) >> (11 + i);
                }
                if (k.n >= 8)
                    ret ^= k.pr_hash_value;
                for (auto i = 0u; i < k.degree.size(); ++i) {
                    ret ^= (hash<int>()(k.degree[i]) << (19 - i)) >> (13 + i);
                }
                return ret;
            default:
                return std::numeric_limits<size_t>::max();
        }
    }
};
