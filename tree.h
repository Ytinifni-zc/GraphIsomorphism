//
// Created by inFinity on 2019-09-17.
//

#ifndef GI_TREE_H
#define GI_TREE_H

#include <vector>
#include <queue>
#include <algorithm>
#include <unordered_set>
#include <set>
#include <bitset>

#include "Type.h"

using Code = std::vector<bool>;
static auto cmp = [](const Code& l, const Code& r) {
    return l.size() < r.size();
};

inline int ij2idx(UInt32 n, int i, int j) {
    if (i > j) std::swap(i, j);
    return n * i + j - i * (i + 1) / 2;
}

inline int get(const std::vector<int>& adj, UInt32 n, int i, int j) {
    return adj[ij2idx(n, i, j)];
}

class Tree {
    using Ngh = std::unordered_set<UInt32>;
private:
    void setCenter() {
        using ChildParent = std::pair<UInt32, UInt32>;
        std::queue<ChildParent> q;

        std::vector<Ngh> tmp_nghs = nghs;

        for (auto i = 0u; i < n; ++i) {
            if (tmp_nghs[i].size() == 1) {
                q.emplace(i, *(tmp_nghs[i].begin()));
            }
//            std::cout << tmp_nghs[i].size() << std::endl;
        }
        center = q.back().first;
        while (!q.empty()) {
            const auto& [c, p] = q.front();
            q.pop();
            tmp_nghs[p].erase(c);
            if (tmp_nghs[p].size() == 1) {
                q.emplace(p, *(tmp_nghs[p].begin()));
            }
            if (!q.empty())
                center = q.back().first;
        }
        getChildren();
    }

    void getChildren(UInt32 r) {
        for (auto& c: nghs[r]) {
            nghs[c].erase(r);
            getChildren(c);
        }
    }

    void getChildren() {
        getChildren(center);
    }

    Code encode(UInt32 r) const {
        if (nghs[r].empty())
            return Code{false, true};
        std::vector<Code> child_codes;
        for (const auto& c: nghs[r]) {
            child_codes.emplace_back(encode(c));
        }
        std::sort(child_codes.begin(), child_codes.end(), cmp);
        Code code{false};
        for (auto &cd: child_codes) {
            code.insert(code.end(), cd.begin(), cd.end());
        }
        code.emplace_back(true);
        return code;
    }

public:
    UInt32 n=0;
    UInt32 center = 0;
    std::vector<Ngh> nghs;

    Tree() = default;

    Tree(UInt32 n, const std::vector<int>& adj) : n(n) {
        nghs.resize(n);
        for (auto i = 0; i < n - 1; ++i) {
            for (auto j = i + 1; j < n; ++j) {
                if (get(adj, n, i, j)) {
                    nghs[i].insert(j);
                    nghs[j].insert(i);
                }
            }
        }
        setCenter();
    }

    auto encode() const {
        return encode(center);
    }

    auto hash() {
        std::bitset<64> bit_code;
        auto code = encode();
        for (auto i = 0u; i < code.size(); ++i) {
            bit_code[i] = code[i];
        }
        return std::hash<std::bitset<64>>()(bit_code);
    }

    void printTree(std::ostream& os, UInt32 r) const {
        os << ", " << r;
        if (nghs[r].empty()) return;
        os << ": [";
        for (auto &c: nghs[r])
            printTree(os, c);
        os << "]";
    }

    friend std::ostream &operator<<(std::ostream &os, const Tree &t) {
        t.printTree(os, t.center);
        os << "{";
        auto code = t.encode();
        for (auto cd: code) {
            if (cd)
                os << 1;
            else
                os << 0;
        }
        os << "}";
        return os;
    }

};

#endif //GI_TREE_H
