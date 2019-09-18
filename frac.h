//
// Created by inFinity on 2019-08-19.
//

#ifndef GI_FRAC_H
#define GI_FRAC_H

#include <iostream>

using namespace std;
using int64 = long long;

class Fraction {
private:
    // Calculates the greates common divisor with
    // Euclid's algorithm
    // both arguments have to be positive
    int64 gcd(int64 a, int64 b) {
        while (a != b) {
            if (a > b) {
                a -= b;
            } else {
                b -= a;
            }
        }
        return a;
    }

public:
    int64 numerator, denominator;

    Fraction() {
        numerator = 0;
        denominator = 1;
    }

    Fraction(int64 n, int64 d) {
        if (d == 0) {
            cerr << "Denominator may not be 0." << endl;
            exit(0);
        } else if (n == 0) {
            numerator = 0;
            denominator = 1;
        } else {
            int sign = 1;
            if (n < 0) {
                sign *= -1;
                n *= -1;
            }
            if (d < 0) {
                sign *= -1;
                d *= -1;
            }

            int64 tmp = gcd(n, d);
            numerator = n / tmp * sign;
            denominator = d / tmp;
        }
    }

    operator int() { return (numerator) / denominator; }

    operator float() { return ((float) numerator) / denominator; }

    operator double() { return ((double) numerator) / denominator; }

    size_t hash() const {
        using std::hash;
        return (hash<int64>()(numerator) << 7) ^ (hash<int64>()(denominator) >> 3);
    }
};

Fraction operator+(const Fraction &lhs, const Fraction &rhs) {
    Fraction tmp(lhs.numerator * rhs.denominator
                 + rhs.numerator * lhs.denominator,
                 lhs.denominator * rhs.denominator);
    return tmp;
}

Fraction operator+(const Fraction &lhs, int rhs) {
    Fraction tmp(rhs, 1);
    return lhs + tmp;
}

Fraction operator+=(Fraction &lhs, const Fraction &rhs) {
    Fraction tmp(lhs.numerator * rhs.denominator
                 + rhs.numerator * lhs.denominator,
                 lhs.denominator * rhs.denominator);
    lhs = tmp;
    return lhs;
}

Fraction operator+=(const Fraction &lhs, int rhs) {
    Fraction tmp(rhs, 1);
    lhs += tmp;
    return lhs;
}

Fraction operator-(const Fraction &lhs, const Fraction &rhs) {
    Fraction tmp(lhs.numerator * rhs.denominator
                 - rhs.numerator * lhs.denominator,
                 lhs.denominator * rhs.denominator);
    return tmp;
}

Fraction operator-(const Fraction &lhs, int rhs) {
    Fraction tmp(rhs, 1);
    return lhs - tmp;
}

Fraction operator-=(Fraction &lhs, const Fraction &rhs) {
    Fraction tmp(lhs.numerator * rhs.denominator
                 - rhs.numerator * lhs.denominator,
                 lhs.denominator * rhs.denominator);
    lhs = tmp;
    return lhs;
}

Fraction operator-=(const Fraction &lhs, int rhs) {
    Fraction tmp(rhs, 1);
    lhs -= tmp;
    return lhs;
}

Fraction operator*(const Fraction &lhs, const Fraction &rhs) {
    Fraction tmp(lhs.numerator * rhs.numerator,
                 lhs.denominator * rhs.denominator);
    return tmp;
}

Fraction operator*=(Fraction &lhs, const Fraction &rhs) {
    Fraction tmp(lhs.numerator * rhs.numerator,
                 lhs.denominator * rhs.denominator);
    lhs = tmp;
    return lhs;
}

Fraction operator*(int lhs, const Fraction &rhs) {
    Fraction tmp(lhs * rhs.numerator, rhs.denominator);
    return tmp;
}

Fraction operator*(const Fraction &rhs, int lhs) {
    Fraction tmp(lhs * rhs.numerator, rhs.denominator);
    return tmp;
}

Fraction operator/(const Fraction &lhs, const Fraction &rhs) {
    Fraction tmp(lhs.numerator * rhs.denominator,
                 lhs.denominator * rhs.numerator);
    return tmp;
}

Fraction operator/(const Fraction &lhs, int rhs) {
    assert(rhs != 0);
    Fraction tmp(lhs.numerator,
                 lhs.denominator * rhs);
    return tmp;
}

Fraction operator/=(Fraction &lhs, const Fraction &rhs) {
    Fraction tmp(lhs.numerator * rhs.denominator,
                 lhs.denominator * rhs.numerator);
    lhs = tmp;
    return lhs;
}

Fraction operator/=(Fraction &lhs, int rhs) {
    assert(rhs != 0);
    Fraction tmp(lhs.numerator,
                 lhs.denominator * rhs);
    lhs = tmp;
    return lhs;
}

bool operator<(const Fraction &lhs, const Fraction &rhs) {
    auto tmp = lhs-rhs;
    return tmp.numerator < 0;
}

bool operator>(const Fraction &lhs, const Fraction &rhs) {
    auto tmp = lhs-rhs;
    return tmp.numerator > 0;
}

std::ostream &operator<<(std::ostream &strm, const Fraction &a) {
    if (a.denominator == 1) {
        strm << a.numerator;
    } else {
        strm << a.numerator << "/" << a.denominator;
    }
    return strm;
}

bool operator==(const Fraction &lhs, const Fraction &rhs) {
    return lhs.numerator == rhs.numerator && lhs.denominator == rhs.denominator;
}

bool operator!=(const Fraction &lhs, const Fraction &rhs) {
    return lhs.numerator != rhs.numerator || lhs.denominator != rhs.denominator;
}

#endif //GI_FRAC_H
