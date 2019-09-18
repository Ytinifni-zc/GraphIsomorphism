//
// Created by inFinity on 2019-09-09.
//

#ifndef GI_BLOOMFILTER_H
#define GI_BLOOMFILTER_H

#include <vector>
#include <array>


//basic structure of a bloom filter object
struct BloomFilter {
private:
    uint8_t m_numHashes;
    std::vector<bool> m_bits;

public:
    BloomFilter(size_t size, uint8_t numHashes) : m_bits(size), m_numHashes(numHashes) {}

    void add(const uint8_t *data, std::size_t len) {
        auto hashValues = myhash(data, len);
        for (int n = 0; n < m_numHashes; n++) {
            m_bits[nthHash(n, hashValues[0], hashValues[1], m_bits.size())] = true;
        }
    }

    bool possiblyContains(const uint8_t *data, std::size_t len) const {
        auto hashValues = myhash(data, len);
        for (int n = 0; n < m_numHashes; n++) {
            if (!m_bits[nthHash(n, hashValues[0], hashValues[1], m_bits.size())]) {
                return false;
            }
        }
        return true;
    }

public:
    inline static size_t nthHash(int n, uint64_t hashA, uint64_t hashB, size_t filterSize) {
        return (hashA + n * hashB) % filterSize; // <- not sure if that is OK, perhaps it is.
    }

};


#include <functional>
#include <iostream>
#include <assert.h>

using namespace std;

int main() {
    BloomFilter bf(1024 * 1024, 5);
    std::string s = "12345";
    bf.add((uint8_t *) s.c_str(), s.size());
    bool possiblyContains = bf.possiblyContains((uint8_t *) s.c_str(), s.size());
    cout << "possible: " << possiblyContains << endl;
    assert(possiblyContains);
}

#endif //GI_BLOOMFILTER_H
