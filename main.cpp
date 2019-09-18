#include <iostream>
#include <cstdint>
#include <unordered_map>
#include <unistd.h>
#include <utils/time_cost.hpp>

#include "Matrix.h"
#include "Gang.h"

size_t threadNum() {
    return sysconf(_SC_NPROCESSORS_ONLN);
}

auto ij2idx(int vn, int i, int j) {
    if (i > j) std::swap(i, j);
    return (vn - 1) * i - i * (i + 1) / 2 + j - 1;
}

Matrix int2Matrix(unsigned long adj, int vn) {
    auto bit = [adj](auto idx) {
        return (adj >> idx) & 1;
    };
    Matrix A(vn);
    for (auto i = 0; i < vn; ++i) {
        for (auto j = i + 1; j < vn; ++j) {
            auto idx = ij2idx(vn, i, j);
            if (bit(idx)) {
                A(i, j) = 1;
            }
        }
    }
    return std::move(A);
}

auto hashGraphs(int vn, Matrix::HashFunction hf) {
    auto upper_bound = (1ul << (vn * (vn - 1) / 2)) - 1;
    std::cout << "Upper Bound: " << upper_bound << std::endl;

    using MatMap = std::unordered_map<Matrix, int, MatrixHasher>;
    const auto thread_num = threadNum();
    std::vector<MatMap> maps(thread_num);
    std::vector<unsigned long> connected_nums(thread_num);

    GangUtil gu;
    gu.submit(upper_bound, [&](UInt64 j) {
        auto i = j + 1;
        auto mt = int2Matrix(i, vn);
        if (mt.isConnected()) {
            auto tid = Gang::thread_id;
            connected_nums[tid]++;
            mt.setHash(hf);
            mt.hash();
            maps[tid][mt]++;
        }
    });
//#pragma omp parallel for
//    for (auto i = 1ul; i <= upper_bound; ++i) {
//        auto mt = int2Matrix(i, vn);
//        if (!mt.isConnected()) continue;
//        auto tid = omp_get_thread_num();
//        connected_nums[tid]++;
//        mt.setHash(hf);
//        mt.hash();
//        maps[tid][mt]++;
//    }

    auto reduce = [&](auto maps) {
        auto &m_map = maps.front();
        for (auto i = 1u; i < maps.size(); ++i) {
            for (auto &[k, v]: maps[i]) {
                m_map[k] += v;
            }
        }
        auto connected_num = 0ul;
        for (auto cn: connected_nums)
            connected_num += cn;

        std::cout << "Number: " << m_map.size() << std::endl;
        std::cout << "Total subgraphs: " << connected_num << std::endl;
        for (auto &[k, v]: m_map) {
            std::cout << k << ": " << v << std::endl;
        }
        return m_map;
    };
    return reduce(maps);
}

int main(int argc, char *argv[]) {
    auto vn = 4;
    Matrix::HashFunction hf = Matrix::Eigen;
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
        } else {
            auto eigen_map = hashGraphs(vn, Matrix::Eigen);
            auto bliss_map = hashGraphs(vn, Matrix::Bliss);
            auto evc_map = hashGraphs(vn, Matrix::EVC);
            auto pr_map = hashGraphs(vn, Matrix::PR);
            auto eig_pr_map = hashGraphs(vn, Matrix::PR);
            std::cout << "-----------------------------" << std::endl;
            auto print_value = [&](auto &emap, auto &k) {
                auto flag = false;
                for (auto &[ek, ev]: emap) {
                    if (ek.adj == k.adj) {
                        std::cout << "\t" << ev;
                        flag = true;
                    }
                }
                if (!flag)
                    std::cout << "\t" << 0;
            };
            for (auto &[k, v]: bliss_map) {
                std::cout << k;
                print_value(eigen_map, k);
                std::cout << "\t" << v;
                print_value(evc_map, k);
                print_value(pr_map, k);
                print_value(eig_pr_map, k);
                std::cout << std::endl;
            }
            return 0;
        }
    }
    common_utils::cost([&]() { hashGraphs(vn, hf); });
    return 0;
}