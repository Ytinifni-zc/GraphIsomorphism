//
// Created by inFinity on 2019-09-17.
//

#include <iostream>
#include <cstdlib>
#include <utils/time_cost.hpp>

#include "tree.h"
#include "Matrix.h"
#include "Gang.h"
#include "StopWatch.h"

auto randomTree(Matrix &m, int vn) {
    // Edge num of tree = vn-1
    std::set<int> candidate;
    std::set<int> pre_tree;
    for (auto i = 0u; i < vn; ++i)
        candidate.insert(i);
    auto u = std::rand() % vn;
    candidate.erase(u);
    pre_tree.insert(u);

    for (auto i = 0u; i < vn - 1; ++i) {
        auto can_size = candidate.size();
        auto idx_u = std::rand() % can_size;
        auto it_u = candidate.begin();
        std::advance(it_u, idx_u);
        auto u = *it_u;
        candidate.erase(u);

        auto pre_tree_size = pre_tree.size();
        auto idx_v = std::rand() % pre_tree_size;
        auto it_v = pre_tree.begin();
        std::advance(it_v, idx_v);
        auto v = *it_v;
        pre_tree.insert(u);
#ifdef PRINT_DEBUG
        std::cout << u << " " << v << "\n";
#endif
        m(u, v) = 1;
    }
#ifdef PRINT_DEBUG
    std::cout << "--------------- \n";
#endif
}

auto genRandTrees(int vn, int rand_num) {
    using Matrices = std::vector<Matrix>;
    Matrices mtxes;
    GangUtil gu;
    common_utils::cost([&]() {
        mtxes.resize(rand_num, Matrix(vn));
        gu.submit(rand_num, [&](uint i) {
            auto &mt = mtxes[i];
            randomTree(mt, vn);
        });
    });
    return std::move(mtxes);
}

auto hashGraphs(int vn, Matrix::HashFunction hf, const std::vector<Matrix> &graph_set) {
    using MatMap = std::unordered_map<Matrix, int, MatrixHasher>;

    GangUtil gu;
    const auto thread_num = gu.thread_num;
    std::vector<MatMap> maps(thread_num);

    common_utils::cost([&]() {
        gu.submit(graph_set.size(), [&](UInt64 i) {
            Matrix mt(vn);
            mt.clone(graph_set[i]);
            auto tid = Gang::thread_id;
            mt.setHash(hf);
            mt.hash();
#ifdef PRINT_DEBUG
            std::cout << mt << std::endl;
#endif
            maps[tid][mt]++;
        });
        std::cout << "Map: ";
    });

    auto reduce = [&](auto maps) {
        auto &m_map = maps.front();
        common_utils::cost([&]() {

            for (auto i = 1u; i < maps.size(); ++i) {
                for (auto &[k, v]: maps[i]) {
                    m_map[k] += v;
                }
            }

            std::cout << "Number: " << m_map.size() << std::endl;
            std::cout << "Reduce: ";
        });
        return std::move(m_map);
    };

    return std::move(reduce(maps));
}

int main(int argc, char *argv[]) {
    auto vn = 4;
    auto rand_num = 10000000;
    if (argc > 1) {
        vn = atoi(argv[1]);
        if (argc > 2) {
            rand_num = atoi(argv[2]);
        }
    }
    auto graph_set = genRandTrees(vn, rand_num);
    std::cout << "graph set generated." << std::endl;
    using MatMap = std::unordered_map<Matrix, int, MatrixHasher>;
    MatMap tree_map, bliss_map;
    common_utils::cost([&]() {
        tree_map = hashGraphs(vn, Matrix::TreeHash, graph_set);
        std::cout << "Tree: ";
    });
    common_utils::cost([&]() {
        bliss_map = hashGraphs(vn, Matrix::Bliss, graph_set);
        std::cout << "Bliss: ";
    });

    auto error_rate = [&](const MatMap &map) {
        auto sum = 0;
        double ret = 0;
        for (const auto &it : map) {
            Matrix tmp = it.first;
            tmp.setHash(Matrix::Bliss);
            tmp.blissing();
            auto v = it.second;
            sum += v;
            auto correct = bliss_map[tmp];
            if (v != correct) {
                ret += abs(v - correct);
            }
        }
        std::cout << ret << " in " << sum << std::endl;
        ret /= sum;
        return ret;
    };

    std::cout << "-----------------------------" << std::endl;
    common_utils::cost([&]() {
        std::cout << "Eigen error rate: " << error_rate(tree_map) << std::endl;
    });

    return 0;
}
