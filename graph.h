#pragma once

#include <vector>
#include <cstdint>
#include <fstream>
#include <parallel/algorithm>
#include <iostream>

using uint = uint32_t;
using ulong = uint64_t;

class Graph {
public:
    uint vertex_num;
    ulong edge_num;

    std::vector<ulong> out_offsets;
    std::vector<uint> out_list;
    std::vector<ulong> in_offsets;
    std::vector<uint> in_list;

public:
    explicit Graph(uint n, ulong m) : vertex_num(n), edge_num(m) {
        out_offsets.resize(n + 1, 0ul);
        out_list.resize(m);
        in_offsets.resize(n + 1, 0ul);
        in_list.resize(m);
    }

    void readFromTSV(const char *filename) {

        struct Edge {
            uint s;
            uint d;
        };
        auto outCmp = [](const Edge &lhs, const Edge &rhs) {
            return (lhs.s < rhs.s) || (lhs.s == rhs.s && lhs.d < rhs.d);
        };
        auto inCmp = [](const Edge &lhs, const Edge &rhs) {
            return (lhs.d < rhs.d) || (lhs.d == rhs.d && lhs.s < rhs.s);
        };

        std::vector<Edge> edge_list(edge_num);
        std::ifstream ifs(filename);
        ulong cnt = 0ul;
        while (!ifs.eof()) {
            int s, d;
            ifs >> s >> d;
            if (s < 0) s = -s;
            if (d < 0) d = -d;
            edge_list[cnt++] = {static_cast<uint>(s), static_cast<uint>(d)};
        }
        ifs.close();
        std::cout << "Edge_list loaded." << std::endl;

        auto writeList = [&](std::vector<ulong> &offsets, std::vector<uint> &list, bool in_or_out) {

            if (in_or_out)
                __gnu_parallel::sort(edge_list.begin(), edge_list.end(), inCmp);
            else
                __gnu_parallel::sort(edge_list.begin(), edge_list.end(), outCmp);
            std::cout << "Edge_list sorted." << std::endl;
#pragma omp parallel for
            for (auto i = 0u; i < edge_num; ++i) {
                if (in_or_out)
                    list[i] = edge_list[i].s;
                else
                    list[i] = edge_list[i].d;
                if (i < edge_num - 1) {
                    if (in_or_out)
                        for (auto cv = edge_list[i].d; cv < edge_list[i + 1].d; ++cv) {
                            offsets[cv + 1] = i + 1;
                        }
                    else
                        for (auto cv = edge_list[i].s; cv < edge_list[i + 1].s; ++cv) {
                            offsets[cv + 1] = i + 1;
                        }
                }
            }
            for (auto i = offsets.size() - 1; offsets[i] == 0ul; --i) {
                offsets[i] = edge_num;
            }
        };
        writeList(in_offsets, in_list, true);
        writeList(out_offsets, out_list, false);
    }

    inline auto outDegree(uint v) { return out_offsets[v + 1] - out_offsets[v]; }

    inline auto inDegree(uint v) { return in_offsets[v + 1] - in_offsets[v]; }

    struct Ngh {
        std::vector<uint>::iterator begin;
        ulong len;
    };

    inline auto outList(uint v) {
        return Ngh{out_list.begin() + out_offsets[v], outDegree(v)};
    }

    inline auto inList(uint v) {
        return Ngh{in_list.begin() + in_offsets[v], inDegree(v)};
    }
};
