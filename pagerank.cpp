#include <iostream>
#include <cassert>
#include <cstdlib>
#include <cstdio>
#include <memory>
#include <utils/time_cost.hpp>

#include "graph.h"
#include "Gang.h"

int main(int argc, char *argv[]) {
    assert(argc == 3);
    uint n = atoi(argv[1]);
    char cmd[1024];
    sprintf(cmd, "wc -l %s | awk '{print $1}'", argv[2]);
    char buf[1024];
    auto fp = popen(cmd, "r");
    fgets(buf, sizeof(buf), fp);
    ulong m = atol(buf);
    pclose(fp);
    std::unique_ptr<Graph> g;
    common_utils::cost([&]() {
        std::cout << "Load Graph ( " << n << ", " << m << " )\n";
        g = std::make_unique<Graph>(n, m);
        g->readFromTSV(argv[2]);
        std::cout << "Load Graph from " << argv[2] << ": ";
    });

    for (int i = 0; i < 10; ++i) {
        std::cout << i << ": " << g->outDegree(i) << ": " << g->inDegree(i) << std::endl;
    }

    auto damping = 0.85;
    GangUtil gu;
    common_utils::cost([&]() {
        std::vector<float> iod(g->vertex_num);
        std::vector<float> pr(g->vertex_num);
        gu.submit(g->vertex_num, [&](UInt32 i){
            iod[i] = (g->outDegree(i) == 0) ? 1.0 : 1.0 / g->outDegree(i);
            pr[i] = 0.25;
        });
        for (auto it = 0; it < 5; ++it) {
            gu.submit(g->vertex_num, [&](UInt32 i){
                auto &p = pr[i];
                auto i_ngh = g->inList(i);
                for (auto j = 0; j < i_ngh.len; ++j) {
                    auto u = *(i_ngh.begin + j);
                    p += pr[u] * iod[u];
                }
                p *= damping;
                p += (1.0 - damping) / g->vertex_num;
            });
            std::cout << "Iteration " << it << "." << std::endl;
        }
        using PrV = std::pair<float, uint>;
        std::vector<PrV> prs(g->vertex_num);
        gu.submit(g->vertex_num, [&](UInt32 i) {
            prs[i] = {pr[i], i};
        });
        __gnu_parallel::sort(prs.begin(), prs.end(), [](const PrV &lhs, const PrV &rhs) {
            return lhs.first > rhs.first;
        });

        for (int i = 0; i < 10; ++i) {
            auto &p = prs[i];
            std::cout << p.second << ": " << p.first << std::endl;
        }

        std::cout << "PageRank Finished in ";

    });

    return 0;
}
