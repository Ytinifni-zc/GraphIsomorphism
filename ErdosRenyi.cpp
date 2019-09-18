//
// Created by inFinity on 2019-08-16.
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
//#include "tree.h"

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
//#pragma omp parallel for
//        for (auto i = 0ul; i < rand_num; ++i) {
//            auto &mt = mtxes[i];
//            randomGraph(mt, threshold);
//            mt.isConnected();
//        }
    });
    return std::move(mtxes);
}

auto hashGraphs(int vn, Matrix::HashFunction hf, const std::vector<Matrix> &graph_set) {
    auto upper_bound = (1ul << (vn * (vn - 1) / 2)) - 1;
    std::cout << "Upper Bound: " << upper_bound << std::endl;

    using MatMap = std::unordered_map<Matrix, int, MatrixHasher>;
    const auto thread_num = threadNum();
    std::vector<MatMap> maps(thread_num);
    std::vector<unsigned long> connected_nums(thread_num);
    GangUtil gu;

    common_utils::cost([&]() {
        gu.submit(graph_set.size(), [&](UInt64 i) {
            Matrix mt(vn);
            mt.clone(graph_set[i]);
            if (mt.is_connected) {
                auto tid = Gang::thread_id;
                connected_nums[tid]++;
                mt.setHash(hf);
                mt.hash();
                maps[tid][mt]++;
            }
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
            auto connected_num = std::accumulate(connected_nums.begin(), connected_nums.end(), 0ul);

            std::cout << "Number: " << m_map.size() << std::endl;
            std::cout << "Total subgraphs: " << connected_num << std::endl;
            std::cout << "Reduce: ";
        });
        return std::move(m_map);
    };

    return std::move(reduce(maps));
}

int main(int argc, char *argv[]) {
    auto vn = 4;
    Matrix::HashFunction hf = Matrix::Eigen;
    auto threshold = 0.3;
    auto rand_num = 10000000;
    if (argc > 1) {
        vn = atoi(argv[1]);
        if (argc > 2) {
            auto argv2 = atoi(argv[2]);
            switch (argv2) {
                case 0:
                    hf = Matrix::Eigen;
                    break;
                case 1:
                    hf = Matrix::Bliss;
                    break;
                case 2:
                    hf = Matrix::EVC;
                    break;
                case 3:
                    hf = Matrix::PR;
                    break;
                case 4:
                    hf = Matrix::EigenPR;
                    break;
                default:
                    printf("The 2nd arg should be 0(Eigen), 1(Bliss), 2(EVC), 3(pagerank), 4(Eigen & Pagerank).\n");
                    break;
            }
            if (argc > 3) {
                threshold = atof(argv[3]);
                assert(threshold < 1 && threshold > 0);
                if (argc > 4) {
                    rand_num = atoi(argv[4]);
                }
            }
        }
    }

    auto graph_set = genRandGraphs(vn, threshold, rand_num);
    if (argc == 2)
        common_utils::cost([&]() { hashGraphs(vn, hf, graph_set); });
    else {
        using MatMap = std::unordered_map<Matrix, int, MatrixHasher>;
        MatMap eigen_map, evc_map, pr_map, ep_map, bliss_map;
        common_utils::cost([&]() {
            eigen_map = hashGraphs(vn, Matrix::Eigen, graph_set);
            std::cout << "Eigen: ";
        });
//        common_utils::cost([&]() {
//            evc_map = hashGraphs(vn, Matrix::EVC, graph_set);
//            std::cout << "EVC: ";
//        });
//        common_utils::cost([&]() {
//            pr_map = hashGraphs(vn, Matrix::PR, graph_set);
//            std::cout << "PR: ";
//        });
//        common_utils::cost([&]() {
//            ep_map = hashGraphs(vn, Matrix::EigenPR, graph_set);
//            std::cout << "EigenPR: ";
//        });
        common_utils::cost([&]() {
            bliss_map = hashGraphs(vn, Matrix::Bliss, graph_set);
            std::cout << "Bliss: ";
        });

        auto error_rate = [&](const MatMap &map) {
            auto sum = 0;
            double ret = 0;
            for (const auto &it : map) {
                Matrix tmp = it.first;
//                auto cnt = 0;
//                for (auto e: tmp.adj) {
//                    if (e != 0) cnt+=1;
//                }
//                if (cnt >= vn) continue;
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

//        auto print_value = [&](auto &emap, auto &k) {
//            auto flag = false;
//            for (auto &[ek, ev]: emap) {
//                auto tmp = ek;
//                tmp.blissing();
//                if (tmp.hash_value == k.hash_value) {
//                    std::cout << "\t" << ev;
//                    flag = true;
//                }
//            }
//            if (!flag)
//                std::cout << "\t" << 0;
//        };
        std::cout << "-----------------------------" << std::endl;
        common_utils::cost([&]() {
            std::cout << "Eigen error rate: " << error_rate(eigen_map) << std::endl;
        });
        /*
        for (auto &[k, v]: bliss_map) {
            auto cnt = 0;
            for (auto e: k.adj) {
                if (e != 0) cnt+=1;
            }
            if (cnt >= vn) continue;
            std::cout << k;
//            print_value(eigen_map, k);
//            std::cout << "\t" << v;
//            print_value(evc_map, k);
//            print_value(pr_map, k);
//            print_value(eig_pr_map, k);
            std::cout << std::endl;
        }
         */
    }
    return 0;
}
