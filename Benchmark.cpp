//
// Created by inFinity on 2019-09-09.
//

#include <iostream>
#include <cstdint>
#include <unordered_map>
#include <unistd.h>
#include <omp.h>
#include <utils/time_cost.hpp>
#include <random>

#include "Matrix.h"
#include "Gang.h"

size_t threadNum() {
    return sysconf(_SC_NPROCESSORS_ONLN);
}

auto randomGraph(int vn, float threshold) {
    Matrix A(vn);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, 99);
    for (auto i = 0; i < vn; ++i) {
        for (auto j = i + 1; j < vn; ++j) {
            auto rand = dis(gen);
            if (rand < threshold * 100)
                A(i, j) = 1;
        }
    }
    return A;
}

auto randomGraph(Matrix &m, float threshold) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, 99);
    for (auto i = 0; i < m.n; ++i) {
        for (auto j = i + 1; j < m.n; ++j) {
            auto rand = dis(gen);
            if (rand < threshold * 100)
                m(i, j) = 1;
        }
    }
}

auto genRandGraphs(int vn, float threshold, int rand_num) {
    using Matrices = std::vector<Matrix>;
    Matrices mtxes;
    GangUtil gu;
    common_utils::cost([&]() {
        mtxes.resize(rand_num, Matrix(vn));
        gu.submit(rand_num, [&](uint i) {
            auto &mt = mtxes[i];
            randomGraph(mt, threshold);
            mt.isConnected();
        });
    }, "s");
    return std::move(mtxes);
}

auto benchGraphs(int vn, std::vector<Matrix> &graph_set) {
    auto upper_bound = (1ul << (vn * (vn - 1) / 2)) - 1;
    std::cout << "Upper Bound: " << upper_bound << std::endl;

    const auto thread_num = threadNum();
    using MatMap = std::unordered_map<Matrix, int, MatrixHasher>;
    std::vector<MatMap> eigen_maps(thread_num);
    std::vector<MatMap> bliss_maps(thread_num);
    std::vector<MatMap> bliss_hash_maps(thread_num);
    using MatHashMap = std::unordered_map<size_t, int>;
    std::vector<MatHashMap> eigen_hmaps(thread_num);
    std::vector<MatHashMap> bliss_hmaps(thread_num);


    std::vector<bool> connected_marks(graph_set.size());
    GangUtil gu;
    std::cout << "Size of MatMap: " << sizeof(MatMap) << std::endl;

//    common_utils::cost([&]() {
//        const size_t PROFIT_SIZE = 100000ul;
////        const size_t PROFIT_SIZE = graph_set.size();
//        gu.submit(thread_num, [&](UInt32 i) {
//            auto &m1 = eigen_hmaps[i];
//            auto &m2 = bliss_hmaps[i];
//            for (size_t j = 0u; j < PROFIT_SIZE; ++j) {
//                m1.emplace(j, 1);
//                m2.emplace(j, 1);
//            }
//            m1.clear();
//            m2.clear();
//        }, 1);
//        std::cout << "Hashmap profit: ";
//    }, "s");

    common_utils::cost([&]() {
        common_utils::cost([&]() {
            gu.submit(graph_set.size(), [&](UInt64 i) {
                Matrix &mt = graph_set[i];
                if (mt.is_connected) {
                    connected_marks[i] = true;
                }
            });
            std::cout << "Make Connecting Filter: ";
        }, "s");

        common_utils::cost([&]() {
            gu.submit(graph_set.size(), [&](UInt64 i) {
                Matrix &mt = graph_set[i];
                if (connected_marks[i]) {
                    mt.setHash(Matrix::Eigen);
                    mt.hash();
                    auto &map = eigen_maps[Gang::thread_id];
                    map[mt]++;
//                    auto &map = eigen_hmaps[Gang::thread_id];
//                    map[mt.hash_value]++;
                }
            });
            std::cout << "Hash Eigen: ";
        }, "s");

        common_utils::cost([&]() {
            gu.submit(graph_set.size(), [&](UInt64 i) {
                Matrix &mt = graph_set[i];
                if (connected_marks[i]) {
                    mt.setHash(Matrix::BlissHash);
                    mt.hash();
                    auto &map = bliss_hash_maps[Gang::thread_id];
                    map[mt]++;
//                    auto &map = bliss_hmaps[Gang::thread_id];
//                    map[mt.hash_value]++;
                }
            });
            std::cout << "Hash BlissHash: ";
        }, "s");

        common_utils::cost([&]() {
            gu.submit(graph_set.size(), [&](UInt64 i) {
                Matrix &mt = graph_set[i];
                if (connected_marks[i]) {
                    mt.setHash(Matrix::Bliss);
                    mt.hash();
                    auto &map = bliss_maps[Gang::thread_id];
                    map[mt]++;
//                    auto &map = bliss_hmaps[Gang::thread_id];
//                    map[mt.hash_value]++;
                }
            });
            std::cout << "Hash Bliss: ";
        }, "s");

        std::cout << "Map: ";
    }, "s");

    auto reduce = [&](auto &maps) {
        auto &m_map = maps.front();
        common_utils::cost([&]() {

            for (auto i = 1u; i < maps.size(); ++i) {
                for (auto &[k, v]: maps[i]) {
                    m_map[k] += v;
                }
            }

            std::cout << "Number: " << m_map.size() << " ";
        }, "s");
//        return std::move(m_map);
    };


    common_utils::cost([&]() {
        reduce(eigen_maps);
//        reduce(eigen_hmaps);
        std::cout << "Reduce EigenMap: ";
    }, "s");
    common_utils::cost([&]() {
        reduce(bliss_hash_maps);
//        reduce(bliss_hmaps);
        std::cout << "Reduce BlissHashMap: ";
    }, "s");

    common_utils::cost([&]() {
        reduce(bliss_maps);
//        reduce(bliss_hmaps);
        std::cout << "Reduce BlissMap: ";
    }, "s");

}

int main(int argc, char *argv[]) {
    auto vn = 4;
    auto threshold = 0.3;
    auto rand_num = 10000000;
    if (argc > 1) {
        vn = atoi(argv[1]);
        if (argc > 2) {
            threshold = atof(argv[2]);
            assert(threshold < 1 && threshold > 0);
            if (argc > 3) {
                rand_num = atoi(argv[3]);
            }
        }
    }

    auto graph_set = genRandGraphs(vn, threshold, rand_num);
    if (argc != 2) {
        benchGraphs(vn, graph_set);
    }
    return 0;
}
