#include <cfloat>
#include <cmath>
#include <iostream>
#include <limits>
#include <assert.h>
#include <cstring>
#include <vector>

using namespace std;

unsigned f_as_uint(float f) {
    unsigned u;  memcpy(&u, &f, sizeof(unsigned));
    return u;
}

template <typename T>
T macheps() {
    T eps = 1;
    while (1 + eps != 1) {
        eps /= 2;
        if (1 + eps / 2 == 1)
            break;
    }
    return eps;
}

template <typename T>
int mant_dig() {
    T a = 1;
    T b = 0.5;
    T c = a + b;
    int i = 0;
    while (a != c) {
        i++;
        b = b / 2;
        c = a + b;
    }
    return i + 1;
}

template <typename T>
pair<int, int> exp_dig() {
    T a = 1.0;
    int max_exp = 0;
    while (a != INFINITY) {
        a *= 2;
        max_exp++;
    }
    a = 1.0;
    int min_exp = 0;
    // or we can subtract size of significand from final answer.
    while (a > std::numeric_limits<T>::min()) {
        a /= 2.0;
        min_exp--;
    }
    return pair<int, int>(min_exp + 1, max_exp);
}

void compare(std::pair<string, float>& p1, std::pair<string, float>& p2) {
    if (p1.second < p2.second)
        std::cout << p1.first << " is less than " << p2.first << endl; 
    else if (p1.second > p2.second)
        std::cout << p2.first << " is less than " << p1.first << endl; 
    else if (p1.second == p2.second)
        std::cout << p1.first << " is equal to " << p2.first << endl; 
    else
        std::cout << p1.first << " is not equal to " << p2.first << endl; 
}

int main() {
    float eps32 = macheps<float>();
    double eps64 = macheps<double>();
    
    // 1
    assert(eps32 == numeric_limits<float>::epsilon());
    assert(eps64 == numeric_limits<double>::epsilon());

    // 2
    std::cout << std::numeric_limits<float>::min_exponent << std::endl;

    pair<int, int> fexps = exp_dig<float>();
    pair<int, int> dexps = exp_dig<double>();
    assert(fexps.first == std::numeric_limits<float>::min_exponent &&
            fexps.second == std::numeric_limits<float>::max_exponent);
    assert(dexps.first == std::numeric_limits<double>::min_exponent &&
            dexps.second == std::numeric_limits<double>::max_exponent);

    // 3
    assert(mant_dig<float>() == std::numeric_limits<float>::digits);
    assert(mant_dig<double>() == std::numeric_limits<double>::digits);

    // 4
    auto toCmp = std::vector<std::pair<string, float>> {
        {"1", 1}, {"1 + e/2", 1 + eps32 / 2.0}, {"1 + e", 1 + eps32}, {"1 + e + e/2", 1 + eps32 + eps32/2.0}
    };
    for (size_t i = 0 ; i < toCmp.size(); i++) {
        auto& p1 = toCmp.at(i);
        for (size_t j = i + 1; j < toCmp.size(); j++) {
            auto& p2 = toCmp.at(j);
            compare(p1, p2);
        }
    }

    return 0;
}
