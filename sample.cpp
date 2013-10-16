#include <iostream>
#include <vector>

#include "alglib.hpp"

int main()
{
    std::vector<int> a = { 1, 3, 5, 7, 9 };
    std::vector<int> b = { 3, 9, 5, 6, 7, 8 };
    std::cout << lcs(a,b) << std::endl;
    return 0;
}
