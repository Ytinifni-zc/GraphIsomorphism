//
// Created by inFinity on 2019-08-15.
//

#include <iostream>
#include <cstdint>
#include <unordered_map>
#include <unistd.h>
#include <utils/time_cost.hpp>
#include <fstream>
#include <string>

#include "Matrix.h"
#include "tree.h"

using namespace std;

auto ij2idx(int vn, int i, int j) {
    if (i > j) std::swap(i, j);
    return vn * i - i * (i + 1) / 2 + j;
}

int main(int argc, char* argv[]) {
    Matrix m1(14);
    m1(0, 1) = 1;
    m1(0, 2) = 1;
    m1(1, 3) = 1;
    m1(1, 4) = 1;
    m1(1, 5) = 1;
    m1(2, 6) = 1;
    m1(2, 7) = 1;
    m1(7, 8) = 1;
    m1(8, 9) = 1;
    m1(8, 10) = 1;
    m1(10, 11) = 1;
    m1(10, 12) = 1;
    m1(11, 13) = 1;
    Tree t(m1.n, m1.adj);
    std::cout << t.center << " hash: " << t.hash() << std::endl;

    Matrix m2(14);
    m2(0, 1) = 1;
    m2(0, 2) = 1;
    m2(1, 3) = 1;
    m2(1, 4) = 1;
    m2(2, 5) = 1;
    m2(2, 6) = 1;
    m2(3, 7) = 1;
    m2(5, 8) = 1;
    m2(5, 9) = 1;
    m2(7, 10) = 1;
    m2(7, 11) = 1;
    m2(7, 12) = 1;
    m2(8, 13) = 1;
    Tree t2(m2.n, m2.adj);
    std::cout << t2.center << " hash: " << t2.hash() << std::endl;
    return 0;
}
