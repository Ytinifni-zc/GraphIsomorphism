//
// Created by inFinity on 2019-08-13.
//

#include <iostream>
#include <cstdint>
#include <unordered_map>
#include <unistd.h>
#include <utils/time_cost.hpp>
#include <fstream>
#include <string>

#include "Matrix.h"

auto ij2idx(int vn, int i, int j) {
    if (i > j) std::swap(i, j);
    return vn * i - i * (i + 1) / 2 + j;
}

int main(int argc, char *argv[]) {
    auto hf = Matrix::Eigen;
    if (argc > 1)
        switch (atoi(argv[1])) {
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
                break;
        }
    using string = std::string;
    auto vn = 8;
    string filename = "error.txt";
    if (argc > 2) {
        vn = atoi(argv[2]);
        if (argc > 3)
            filename = string(argv[3]);
    }
    auto gen_mtx = [&](string s) {
        Matrix m(vn, hf);
        for (auto i = 0; i < vn; ++i) {
            for (auto j = i + 1; j < vn; ++j) {
                if (s[ij2idx(vn, i, j)] == '1') {
                    m(i, j) = 1;
                }
            }
        }
//        std::cout << s << std::endl;
//        m.setDegree();
//        m.printMatrix();
        m.hash();
//        std::cout << std::endl;
        return m;
    };
    auto process_string_graph = [&](string s) {
        auto m = gen_mtx(s);
        MatrixHasher mh;
        return mh(m);
    };

    std::ifstream ifs(filename);
    std::unordered_map<Matrix, std::vector<string>, MatrixHasher> map;
    while (!ifs.eof()) {
        string s;
        std::getline(ifs, s);
        if (s.empty()) break;
        auto m = gen_mtx(s);
#ifdef PRINT_DEBUG
        m.printMatrix();
        MatrixHasher mh;
        std::cout << mh(m) << std::endl;
#endif
        map[m].push_back(s);
//        auto hv = process_string_graph(s);
//        map[hv].push_back(s);
    }
    ifs.close();
    std::cout << map.size() << std::endl;
    for (auto&[m, v]: map) {
//        auto &s = v.front();
//        for (auto c: s) {
//            if (c == '1') e++;
//        }
        if (v.size() < 2) continue;
        auto e = 0;
        for (auto d: m.degree) e += d;
        std::cout << e/2 << ": ";
        for (auto &s: v) {
            std::cout << s << " ";
        }
        std::cout << std::endl;
//        std::cout << "[ ";
//        for (auto e: m.eigen) {
//            std::cout << e << " ";
//        }
//        std::cout << "]" << std::endl;
    }
    return 0;
}
