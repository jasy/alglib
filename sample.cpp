#include <iostream>
#include <string>
#include <vector>
#include <functional>
#include <cassert>

#include "alglib.hpp"

int main()
{
    // Greatest Common Divisor
    // Least Common Multiple
    {
        const int a = 2*3*3*7*13*13*13;
        const int b = 2*3*5*7*11*13*13;
        assert(2*3*7*13*13==gcd(a,b));
        assert(2*3*3*5*7*11*13*13*13==lcm(a,b));
    }
    // Modulo Integer
    {
        typedef mint<int,13> mi;
        assert(1==(mi(9)+5).get());
        assert(12==(mi(1)-2).get());
        assert(2==(mi(3)*5).get());
        assert(7==(mi(1)/2).get());
        assert(3==mi(3).pow(4).get());
        assert(6==mi::pow(2,5).get());
        assert(2==mi::c(10,4).get());
        mint<int> a(100000);
        assert(999999937==(a*a).get()); // no overflow
    }
    // Longest Common Subsequence
    assert(3==lcs<std::string>("13579","395678"));
    {
        std::vector<int> a = {  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27 };
        std::vector<int> b = {  2,  3, 30,  4,  5,  6, 31,  7, 32,  9, 10, 11, 12, 33, 14, 15, 16, 34, 17, 19, 20, 22, 35, 24, 36, 37, 25, 26, 27, 38, 39 };
        assert(21==lcs(a,b));
    }
    // Warshall-Floyd
    {
        const int N = 100;
        std::vector<std::vector<int>> a(N,std::vector<int>(N,1000000000));
        for(int i=0; i<N; ++i)
            a[i][i]=0;
        for(int i=0; i<N-1; ++i)
            a[i][i+1]=a[i+1][i]=1;
        wf(a);
        assert(N-1==a[0][N-1]);
        std::vector<std::vector<bool>> b(N,std::vector<bool>(N,false));
        for(int i=0; i<N; ++i)
            b[i][i]=true;
        for(int i=0; i<N-1; ++i)
            b[i][i+1]=b[i+1][i]=true;
        wf<std::logical_and,std::logical_or>(b);
        assert(b[0][N-1]);
    }
    return 0;
}
