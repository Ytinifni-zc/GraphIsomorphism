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
static auto cmp = [](const Code &l, const Code &r) {
    if (l.size() != r.size())
        return l.size() < r.size();
    for (auto i = 0u; i < l.size(); ++i) {
        if (l[i] < r[i]) return true;
        if (l[i] > r[i]) return false;
    }
    return false;
};

inline int ij2idx(UInt32 n, int i, int j) {
    if (i > j) std::swap(i, j);
    return n * i + j - i * (i + 1) / 2;
}

inline int get(const std::vector<int> &adj, UInt32 n, int i, int j) {
    return adj[ij2idx(n, i, j)];
}

class Tree {
    using Ngh = std::unordered_set<UInt32>;
private:
    void setCenter() {
        // Erase
        // https://stackoverflow.com/questions/5055964/centre-node-in-a-tree
        /*
        using ChildParent = std::pair<UInt32, UInt32>;
        using ChildParentQ = std::queue<ChildParent>;
        ChildParentQ q;

        std::vector<Ngh> tmp_nghs = nghs;

        for (auto i = 0u; i < n; ++i) {
            if (tmp_nghs[i].size() == 1) {
                q.emplace(i, *(tmp_nghs[i].begin()));
            }
        }
        bool singleOrDouble = (n - q.size())%2 == 1;

        center = q.back().first;
        auto center2 = q.front().first;
        while (!q.empty()) {
            const auto& [c, p] = q.front();
            q.pop();
            tmp_nghs[p].erase(c);
            if (tmp_nghs[p].size() == 1) {
                q.emplace(p, *(tmp_nghs[p].begin()));
            }
            if (!q.empty()) {
                center = q.back().first;
                if (q.size() == 2)
                    center2 = q.front().first;
            }
        }
        if (!singleOrDouble) {
            // Add Super node
            nghs.push_back({center, center2});
            nghs[center].erase(center2);
            nghs[center2].erase(center);
            center = n;
        }
         */

        auto BFS = [&](UInt32 start) {
            std::vector<bool> visited(n);
            std::queue<UInt32> q;
            q.push(start);
            visited[start] = true;
            auto v2 = q.back();
            while (!q.empty()) {
                auto v = q.front();
                visited[v] = true;
                q.pop();
                for (auto u: nghs[v]) {
                    if (!visited[u])
                        q.push(u);
                }
                if (q.empty())
                    v2 = q.back();
            }
            return v2;
        };
        auto v2 = BFS(0);

        using Log = std::pair<UInt32, UInt32>;
        std::vector<Log> log_stack{{v2, n}};
        std::vector<bool> visited(n);
        std::queue<UInt32> q;
        q.push(v2);
        visited[v2] = true;
        while (!q.empty()) {
            auto v = q.front();
            visited[v] = true;
            q.pop();
            for (auto u: nghs[v]) {
                if (!visited[u]) {
                    q.push(u);
                    log_stack.emplace_back(u, v);
                }
            }
        }
//        std::cout << v2 << std::endl;

        std::vector<UInt32> diameter;
        auto [v3, v3p] = log_stack.back();
        diameter.push_back(v3);
        log_stack.pop_back();
        while (!log_stack.empty()) {
            const auto& [c, p] = log_stack.back();
            if (c == v3p) {
                diameter.push_back(c);
                v3p = p;
            }
            log_stack.pop_back();
        }

        if (diameter.size() % 2)
            center = diameter[diameter.size()/2];
        else {
            auto c1 = diameter[diameter.size()/2];
            auto c2 = diameter[diameter.size()/2-1];
            nghs.push_back({c1, c2});
            nghs[c1].erase(c2);
            nghs[c2].erase(c1);
            center = n;
        }

        getChildren();
    }

    void getChildren(UInt32 r) {
        for (auto &c: nghs[r]) {
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
        for (const auto &c: nghs[r]) {
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
    UInt32 n = 0;
    UInt32 center = 0;

    std::vector<Ngh> nghs;

    Tree() = default;

    Tree(UInt32 n, const std::vector<int> &adj) : n(n) {
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
        return std::hash<std::bitset<64>>()(bit_code) ^ std::hash<UInt32>()(n);
    }

    void printTree(std::ostream &os, UInt32 r) const {
        os << " " << r;
        if (nghs[r].empty()) return;
        os << ":[";
        for (auto &c: nghs[r])
            printTree(os, c);
        os << " ]";
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
