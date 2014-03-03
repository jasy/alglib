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
    // Extended Greatest Common Divisor
    {
        const int a = 2*3*3*7*13*13*13;
        const int b = 2*3*5*7*11*13*13;
        assert(2*3*7*13*13==gcd(a,b));
        assert(2*3*3*5*7*11*13*13*13==lcm(a,b));
        int x=0,y=0;
        assert(2*3*7*13*13==extgcd(a,b,x,y));
        assert(x==24);
        assert(y==-17);
    }
    // Prime Factors
    assert(primes(103)==std::vector<int>({103}));
    assert(primes(2*3*3*7*13*13*13)==std::vector<int>({2,3,3,7,13,13,13}));
    // Integer Square Root
    {
        assert(0==isqrt(0));
        assert(1==isqrt(1));
        assert(1==isqrt(3));
        assert(2==isqrt(4));
        assert(2==isqrt(8));
        assert(3==isqrt(9));
        assert(3999999999ULL==isqrt(15999999999999999999ULL));
        assert(4000000000ULL==isqrt(16000000000000000000ULL));
        assert(4000000000ULL==isqrt(16000000000000000001ULL));
    }
    // Modulo Integer
    {
        typedef mint<int,13> mi;
        assert(1==(mi(9)+5)());
        assert(1==(9+mi(5))());
        assert(12==(mi(1)-2)());
        assert(12==(1-mi(2))());
        assert(2==(mi(3)*5)());
        assert(2==(3*mi(5))());
        assert(7==(mi(1)/2)());
        assert(7==(1/mi(2))());
        assert(3==mi(3).pow(4)());
        assert(6==mi::pow(2,5)());
        assert(2==mi::c(10,4)());
        mint<int> a(100000);
        assert(999999937==(a*a)()); // no overflow
        typedef mint<short,13> ms;
        assert(1==(ms(9)+5)());
        assert(6==ms::pow(2,5)());
        assert(2==ms::c(10,4)());
        mint<long long> b(100000);
        assert(999999937==(b*b)());
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
    // Area of triangle
    assert(0.5==area_of_triangle(1,1,0,0,1,0));
    assert(  6==area_of_triangle(3,4,3,0,0,4));
    // Probability of Complete Gacha
    {
        std::vector<int> a = { 1, 1 };
        assert(3==comp_gacha(a));
        assert(3==comp_gacha_avg(2));
        std::vector<int> b = { 1, 4, 5 };
        auto s145 = comp_gacha(b);
        assert(10.722222 < s145);assert(s145 < 10.722223);
        auto avg50 = comp_gacha_avg(50);
        assert(224.960266 < avg50);assert(avg50 < 224.960267);
    }
    return 0;
}
