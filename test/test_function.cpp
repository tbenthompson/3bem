#include "UnitTest++.h"
#include "function.h"
#include "autocheck/autocheck.hpp"
namespace ac = autocheck;

using namespace tbem;

TEST(Concatenate) {
    ConcatenatedFunction block_fnc = concatenate({
        {1,2},
        {3,4,5},
        {6}
    });
    CHECK_EQUAL(block_fnc.components, 3);
    CHECK_ARRAY_EQUAL(block_fnc.data, (std::vector<double>{1,2,3,4,5,6}), 6);
    CHECK_ARRAY_EQUAL(block_fnc.component_lengths, (std::vector<int>{2,3,1}), 3);
}

TEST(Expand) {
    ConcatenatedFunction block_fnc = concatenate({
        {1,2},
        {6}
    });
    auto expanded = expand(block_fnc);
    CHECK_EQUAL(expanded[0][0], 1);
    CHECK_EQUAL(expanded[0][1], 2);
    CHECK_EQUAL(expanded[1][0], 6);
}

bool expand_concat_identity(const std::vector<std::vector<double>>& A) {
    auto block_fnc = concatenate(A); 
    auto expanded = expand(block_fnc);
    if (A.size() != expanded.size()) {
        return false;
    }
    for (size_t i = 0; i < A.size(); i++) {
        if (A[i].size() != expanded[i].size()) {
            return false;
        }
        for (size_t j = 0; j < A[i].size(); j++) {
            if (fabs(A[i][j] - expanded[i][j]) > 1e-12) {
                return false;
            }
        }
    }
    return true;
}

TEST(ExpandConcatProperty) {
    ac::check<std::vector<std::vector<double>>>(expand_concat_identity);
}

int main(int, char const *[])
{
    return UnitTest::RunAllTests();
}
