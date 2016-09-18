#include <iostream>
#include <string>
#include <vector>
#include <functional>
#include <map>
#include <array>
#include <tuple>
#include <cstdint>
#include <cassert>

#include "alglib.hpp"

#include "gtest/gtest.h"

// Greatest Common Divisor
// Least Common Multiple
// Extended Greatest Common Divisor
TEST(Integer, GCD_LCM_ExGCD)
{
    const int a = 2*3*3*7*13*13*13;
    const int b = 2*3*5*7*11*13*13;
    EXPECT_EQ(2*3*7*13*13, gcd(a,b));
    EXPECT_EQ(2*3*3*5*7*11*13*13*13, lcm(a,b));
    int x=0,y=0;
    EXPECT_EQ(2*3*7*13*13, extgcd(a,b,x,y));
    EXPECT_EQ(x, 24);
    EXPECT_EQ(y, -17);
}

// Prime
TEST(Integer, Prime)
{
    EXPECT_TRUE (prime(2));
    EXPECT_TRUE (prime(3));
    EXPECT_FALSE(prime(4));
    EXPECT_TRUE (prime(5));
    EXPECT_FALSE(prime(6));
    EXPECT_TRUE (prime(7));
    EXPECT_FALSE(prime(8));
    EXPECT_FALSE(prime(9));
    EXPECT_TRUE (prime(9999991));
    EXPECT_FALSE(prime(99999991));
    EXPECT_TRUE (prime(129402307));
    EXPECT_TRUE (prime(282437925089LL));
}

// Prime Factors
TEST(Integer, PrimeFactors)
{
    EXPECT_EQ(std::vector<int>({103}), primes(103));
    EXPECT_EQ(std::vector<int>({2,3,3,7,13,13,13}), primes(2*3*3*7*13*13*13));
}

// Divisors
TEST(Integer, Divisors)
{
    EXPECT_EQ(std::vector<int>({1,129402307}), divisors(129402307));
    EXPECT_EQ(std::vector<int>({1,2,3,5,6,10,15,30}), divisors(2*3*5));
    EXPECT_EQ(std::vector<int>({1,2,3,4,6,9,12,18,36}), divisors(2*2*3*3));
    EXPECT_EQ(std::vector<int>({1,2,3,4,5,6,10,12,15,20,30,60}), divisors(2*2*3*5));
}

// Euler's phi function
TEST(Integer, EulerPhi)
{
    EXPECT_EQ(1, euler_phi(1));
    EXPECT_EQ(1, euler_phi(2));
    EXPECT_EQ(2, euler_phi(3));
    EXPECT_EQ(2, euler_phi(4));
    EXPECT_EQ(4, euler_phi(5));
    EXPECT_EQ(2, euler_phi(6));
    EXPECT_EQ(6, euler_phi(7));
    EXPECT_EQ(4, euler_phi(8));
    EXPECT_EQ(6, euler_phi(9));
    EXPECT_EQ(4, euler_phi(10));
    EXPECT_EQ(10, euler_phi(11));
    EXPECT_EQ(4, euler_phi(12));
    EXPECT_EQ(12, euler_phi(13));
    EXPECT_EQ(6, euler_phi(14));
    EXPECT_EQ(8, euler_phi(15));
    EXPECT_EQ(8, euler_phi(16));
    EXPECT_EQ(2592, euler_phi(5186));
    EXPECT_EQ(2592, euler_phi(5187));
    EXPECT_EQ(2592, euler_phi(5188));
}

// Integer Square Root
TEST(Integer, SquareRoot)
{
    EXPECT_EQ(0, isqrt(0));
    EXPECT_EQ(1, isqrt(1));
    EXPECT_EQ(1, isqrt(3));
    EXPECT_EQ(2, isqrt(4));
    EXPECT_EQ(2, isqrt(8));
    EXPECT_EQ(3, isqrt(9));
    EXPECT_EQ(3999999999ULL, isqrt(15999999999999999999ULL));
    EXPECT_EQ(4000000000ULL, isqrt(16000000000000000000ULL));
    EXPECT_EQ(4000000000ULL, isqrt(16000000000000000001ULL));
}

// Combination
TEST(Integer, Combination)
{
    EXPECT_EQ(1, C(1,0));
    EXPECT_EQ(1, C(1,1));
    EXPECT_EQ((10*9*8*7*6)/(5*4*3*2*1), C(10,5));
    EXPECT_EQ((20*19)/(2*1), C(20,2));
    EXPECT_EQ((20*19*18)/(3*2*1), C(20,17));
    EXPECT_EQ(5200300ULL, C(25ULL,12ULL));
    const int N=25;
    const auto comb = Cs<int>(N); // No overflow
    typedef unsigned long long T;
    for(T n=0; n<=N; ++n) for(T k=0; k<=n; ++k){ EXPECT_EQ(C(n,k),T(comb[n][k])); }
}

// Modulo Integer
TEST(Integer, ModuloInteger)
{
    typedef mint<int,13> mi;
    EXPECT_EQ(1, (mi(9)+5)());
    EXPECT_EQ(1, (9+mi(5))());
    EXPECT_EQ(12, (mi(1)-2)());
    EXPECT_EQ(12, (1-mi(2))());
    EXPECT_EQ(2, (mi(3)*5)());
    EXPECT_EQ(2, (3*mi(5))());
    EXPECT_EQ(7, (mi(1)/2)());
    EXPECT_EQ(7, (1/mi(2))());
    EXPECT_EQ(3, mi(3).pow(4)());
    EXPECT_EQ(6, mi::pow(2,5)());
    EXPECT_EQ(2, mi::c(10,4)());
    mint<int> a(100000);
    EXPECT_EQ(999999937, (a*a)()); // no overflow
    typedef mint<short,13> ms;
    EXPECT_EQ(1, (ms(9)+5)());
    EXPECT_EQ(6, ms::pow(2,5)());
    EXPECT_EQ(2, ms::c(10,4)());
    mint<long long> b(100000);
    EXPECT_EQ(999999937, (b*b)());
}

// Square Matrix
TEST(SquareMatrix, Integer)
{
    enum{ N=3 };
    SquareMatrix<int> m(N);
    for(int i=0; i<N; ++i)
        m(i,i) = 1;
    m*=m*m;
    for(int r=0; r<N; ++r) for(int c=0; c<N; ++c)
    { EXPECT_EQ((r==c?1:0), m(r,c)); }
    {
        auto t = m.pow(1000000000);
        for(int r=0; r<N; ++r) for(int c=0; c<N; ++c)
        { EXPECT_EQ((r==c?1:0), t(r,c)); }
    }
    for(int r=0; r<N; ++r) for(int c=0; c<N; ++c)
        m(r,c) = r+1;
    {
        auto t = m*m;
        for(int r=0; r<N; ++r) for(int c=0; c<N; ++c)
        { EXPECT_EQ(N*(N+1)/2*(r+1), t(r,c)); }
    }
    {
        auto t = m.pow(2);
        for(int r=0; r<N; ++r) for(int c=0; c<N; ++c)
        { EXPECT_EQ(N*(N+1)/2*(r+1), t(r,c)); }
    }
    std::vector<int> v(N);
    for(int i=0; i<N; ++i) v[i]=i+1;
    {
        auto t = m*v;
        for(int i=0; i<N; ++i)
        { EXPECT_EQ(N*(N+1)/2*(i+1), t[i]); }
    }
}

TEST(SquareMatrix, bit)
{
    typedef uint32_t T;
    enum{ N=30 };
    const std::vector<T> A = {11627,5078,8394,6412,10346,3086,3933,668,9879,11739,4501,6108,12336,8771,2768,2438,2153,7047,5476,313,1264,369,12070,10743,10663,747,370,4671,5235,3439};
    const std::vector<T> C = {114,3613,3271,5032,11241,6961,3628,150,12191,2396,7638,3046,11594,8162,11136,786,9878,2356,11660,1070,3649,10882,9746,1415,3307,7077,9319,9981,3437,544};
    const T MulUnit = -1;
    SquareMatrix<T,std::bit_xor,std::bit_and,0,MulUnit> m(N);
    for(int i=0; i<N-1; ++i)
        m(i,i+1)=MulUnit;
    for(int i=0; i<N; ++i)
        m(N-1,i)=C[N-1-i];
    EXPECT_EQ(T(2148), (m.pow(999999999-N)*A)[N-1]);
}

// Fenwick Tree (Binary Indexed Tree)
TEST(Fenwick, standard)
{
    enum{ N=10 };
    fenwick<int> t(N);
    t.add(0,1);
    t.add(N/2,1);
    EXPECT_EQ(0, t.sum(0));
    EXPECT_EQ(1, t.sum(1));
    EXPECT_EQ(1, t.sum(N/2));
    EXPECT_EQ(2, t.sum(N/2+1));
    EXPECT_EQ(2, t.sum(N));
    EXPECT_EQ(size_t(0), t.lower_bound(0));
    EXPECT_EQ(size_t(0), t.lower_bound(1));
    EXPECT_EQ(size_t(N/2), t.lower_bound(2));
    EXPECT_EQ(size_t(N), t.lower_bound(3));
    t.add(N/2,-1);
    t.add(0,-1);
    for(int i=0; i<N; ++i) t.add(i,i+1);
    EXPECT_EQ(N*(N+1)/2, t.sum(N));
    for(int i=0; i<N; ++i){ EXPECT_EQ(i*(i+1)/2, t.sum(i)); }
    for(int i=0; i<N; ++i) t.add(i,N-(i+1));
    EXPECT_EQ(N*N, t.sum(0,N));
    for(int i=0; i<N; ++i){ EXPECT_EQ(N*(N-i), t.sum(i,N)); }
}

TEST(Fenwick, range)
{
    enum{ N=10 };
    fenwick_range<int> t(N);
    for(int i=0; i<N; ++i) t.add(i,N,1);
    EXPECT_EQ(N*(N+1)/2, t.sum(N));
    for(int i=0; i<N; ++i){ EXPECT_EQ(i*(i+1)/2, t.sum(i)); }
    for(int i=0; i<N; ++i) t.add(0,i,1);
    EXPECT_EQ(N*N, t.sum(0,N));
    for(int i=0; i<N; ++i){ EXPECT_EQ(N*(N-i), t.sum(i,N)); }
}

TEST(Fenwick, 2D)
{
    enum{ N=3, M=5 };
    fenwick2D<int> t(N,M);
    for(int i=0; i<N; ++i) for(int j=0; j<M; ++j) t.add(i,j,1);
    EXPECT_EQ(N*M, t.sum(N,M));
    for(int i=0; i<N; ++i) for(int j=0; j<M; ++j){ EXPECT_EQ(i*j, t.sum(i,j)); }
    for(int i1=0; i1<N; ++i1) for(int j1=0; j1<M; ++j1)
        for(int i2=i1; i2<N; ++i2) for(int j2=j1; j2<M; ++j2)
            { EXPECT_EQ((i2-i1)*(j2-j1), t.sum(i1,j1,i2,j2)); }
}

// Levenshtein Distance (Edit Distance)
TEST(Subsequence, EditDistance)
{
    EXPECT_EQ(0, ld<std::string>("",""));
    EXPECT_EQ(0, ld<std::string>("abc","abc"));
    EXPECT_EQ(4, ld<std::string>("abcdef","abfghi"));
    EXPECT_EQ(2, ld<std::string>("abcdef","bcdefg"));
    {
        std::vector<int> a = {  1,  2,  3,  4,  5,  6,  7,  8,  9, 10 };
        std::vector<int> b = {  2,  3, 30,  4,  5,  6, 31,  7, 32,  9, 10 };
        EXPECT_EQ(4, ld(a,b));
    }
}

// Longest Common Subsequence
TEST(Subsequence, LongestCommon)
{
    EXPECT_EQ(3, lcs<std::string>("13579","395678"));
    {
        std::vector<int> a = {  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27 };
        std::vector<int> b = {  2,  3, 30,  4,  5,  6, 31,  7, 32,  9, 10, 11, 12, 33, 14, 15, 16, 34, 17, 19, 20, 22, 35, 24, 36, 37, 25, 26, 27, 38, 39 };
        EXPECT_EQ(21, lcs(a,b));
    }
}

// Longest Increasing Subsequence
TEST(Subsequence, LongestIncreasing)
{
    EXPECT_EQ(4, lis(std::string("135246")));
    EXPECT_EQ(1, lis(std::string("54321")));
    EXPECT_EQ(7, lis(std::string("1234567")));
    EXPECT_EQ(6, lis(std::vector<int>({1,3,4,5,6,7,2})));
}

// Warshall-Floyd
TEST(Graph, WarshallFloyd)
{
    const int N = 10;
    std::vector<std::vector<int>> a(N,std::vector<int>(N,1000000000));
    for(int i=0; i<N; ++i)
        a[i][i]=0;
    for(int i=0; i<N-1; ++i)
        a[i][i+1]=a[i+1][i]=1;
    wf(a);
    EXPECT_EQ(N-1, a[0][N-1]);
    std::vector<std::vector<bool>> b(N,std::vector<bool>(N,false));
    for(int i=0; i<N; ++i)
        b[i][i]=true;
    for(int i=0; i<N-1; ++i)
        b[i][i+1]=b[i+1][i]=true;
    wf<std::logical_and,std::logical_or>(b);
    EXPECT_TRUE(b[0][N-1]);
}

// Shortest Path Faster Algorithm
TEST(Graph, SPFA)
{
    enum{ N=3 };
    typedef int Cost;
    typedef int Vertex;
    typedef SPFA<Cost,Vertex> S;
    std::vector<std::vector<S::Edge>> edges(N);
    for(int i=0; i<N; ++i)
        for(int j=0; j<i; ++j)
            edges[i].emplace_back((i-j)*(i-j),j);
    auto cost = S::spfa(edges, N-1);
    EXPECT_EQ(N, cost.size());
    for(int i=0; i<N; ++i){ EXPECT_EQ(N-1-i, cost[i]); }
}

// Dijkstra's algorithm
TEST(Graph, Dijkstra)
{
    enum{ N=3 };
    typedef int Cost;
    typedef int Vertex;
    typedef Dijkstra<Cost,Vertex> D;
    std::vector<std::vector<D::Edge>> edges(N);
    for(int i=0; i<N; ++i)
        for(int j=0; j<i; ++j)
            edges[i].emplace_back((i-j)*(i-j),j);
    auto cost = D::dijkstra(edges, N-1);
    EXPECT_EQ(N, cost.size());
    for(int i=0; i<N; ++i){ EXPECT_EQ(N-1-i, cost[i]); }
}

// Lowest Common Ancestor
TEST(Graph, LCA)
{
    std::vector<std::pair<int,int>> edges1 = {{2,0},{1,0},{1,3},{1,4},{2,5},{2,6}};
    const LCA lca1(edges1);
    EXPECT_EQ(2, lca1.distance(3,4));
    EXPECT_EQ(2, lca1.distance(0,5));
    EXPECT_EQ(4, lca1.distance(4,5));
    EXPECT_EQ(4, lca1.distance(3,6));
    EXPECT_EQ(3, lca1.distance(4,2));
    std::vector<std::pair<int,int>> edges2 = {{0,1},{0,2},{0,3},{3,4}};
    const LCA lca2(edges2);
    EXPECT_EQ(2, lca2.distance(1,2));
    EXPECT_EQ(2, lca2.distance(1,3));
    EXPECT_EQ(3, lca2.distance(1,4));
}

// Union Find
TEST(Graph, UnionFind)
{
    {
        UnionFind<int> uf(6);
        EXPECT_FALSE(uf.same(2,5));
        std::vector<std::pair<int,int>> edges = {{0,4},{0,5},{1,3},{1,4},{2,3}};
        for(auto& e: edges)
            uf.unite(e.first,e.second);
        for(auto& e: edges)
            EXPECT_TRUE (uf.same(e.first,e.second));
        EXPECT_TRUE (uf.same(2,5));
    }
    {
        UnionFind<int> uf(6);
        EXPECT_FALSE(uf.same(2,5));
        EXPECT_FALSE(uf.same(1,2));
        EXPECT_FALSE(uf.same(4,5));
        std::vector<std::pair<int,int>> edges = {{0,4},{0,5},{1,3},{2,3}};
        for(auto& e: edges)
            uf.unite(e.first,e.second);
        for(auto& e: edges)
            EXPECT_TRUE (uf.same(e.first,e.second));
        EXPECT_FALSE(uf.same(2,5));
        EXPECT_TRUE (uf.same(1,2));
        EXPECT_TRUE (uf.same(4,5));
    }
}

// Kruskal's algorithm
TEST(Graph, Kruskal)
{
    typedef int Cost;
    typedef int Vertex;
    typedef std::tuple<Cost,Vertex,Vertex> Edge;
    std::vector<Edge> edges = {
        Edge{ 40,0,1},
        Edge{ 50,0,2},
        Edge{ 30,0,3},
        Edge{ 70,0,4},
        Edge{ 70,0,5},
        Edge{ 80,0,6},
        Edge{ 80,0,7},
        Edge{ 40,1,2},
        Edge{ 50,1,3},
        Edge{ 60,1,4},
        Edge{ 90,2,5},
        Edge{ 80,3,4},
        Edge{110,4,5},
        Edge{ 60,5,6},
        Edge{ 50,6,7},
    };
    Cost t=0;
    for(auto& e: kruskal(edges, Vertex(8)))
        t += std::get<0>(e);
    EXPECT_EQ(350, t);
}

// Dinic's algorithm
TEST(Graph, Dinic)
{
    typedef int Capacity;
    typedef int Vertex;
    typedef std::tuple<Capacity,Vertex,Vertex> Edge;
    {
        std::vector<Edge> edges = {
            Edge{1,0,1},
            Edge{1,1,2},
            Edge{1,1,3},
            Edge{1,2,4},
            Edge{1,3,4},
        };
        const Vertex N = 5;
        EXPECT_EQ(1, dinic(edges,N,0,N-1));
    }
    {
        std::vector<Edge> edges = {
            Edge{1,0,1},
            Edge{1,0,2},
            Edge{1,1,3},
            Edge{1,2,3},
            Edge{1,3,4},
        };
        const Vertex N = 5;
        EXPECT_EQ(1, dinic(edges,N,0,N-1));
    }
    {
        std::vector<Edge> edges = {
            Edge{1,0,1},
            Edge{1,0,2},
            Edge{1,0,3},
            Edge{1,0,4},
            Edge{1,1,5},
            Edge{1,2,5},
            Edge{1,5,6},
            Edge{1,6,7},
            Edge{1,6,8},
            Edge{1,3,9},
            Edge{1,4,9},
            Edge{1,7,10},
            Edge{1,8,10},
            Edge{1,9,10},
        };
        const Vertex N = 11;
        EXPECT_EQ(2, dinic(edges,N,0,N-1));
    }
}

// Bipartite Matching
TEST(Graph, BipartiteMatching)
{
    std::vector<std::pair<int,int>> edges1 = {{0,3},{0,4},{0,5},{1,3},{1,4},{2,3}};
    EXPECT_EQ(3, bm(edges1,6));
    std::multimap<int,int> edges2 = {{0,4},{0,6},{0,7},{1,4},{1,5},{1,8},{2,5},{2,8},{3,6}};
    EXPECT_EQ(4, bm(edges2,9));
    typedef std::tuple<int,int> T3;
    std::array<T3,1> edges3 = {{T3{0,1}}};
    EXPECT_EQ(1, bm(edges3,2));
}

// Argument of Vector 2D
TEST(Geometry, Argument)
{
    typedef Vec2D<int> T;
    const std::vector<T> p = {T(1,0),T(1,1),T(0,1),T(-1,1),T(-1,0),T(-1,-1),T(0,-1),T(1,-1)};
    const int N = p.size();
    for(int i=0; i<N; ++i){ EXPECT_EQ((i<N/2) ,  p[i].upperArg()); }
    for(int i=0; i<N; ++i)
        for(int j=0; j<N; ++j)
        { EXPECT_EQ((i<j) ,  T::cmpArg(p[i],p[j])); }
    std::vector<T> q(p.rbegin(),p.rend());
    for(int i=0; i<N; ++i){ EXPECT_EQ(p[i], q[N-1-i]); }
    EXPECT_NE(p, q);
    std::sort(q.begin(),q.end(),T::cmpArg);
    EXPECT_EQ(p, q);
}

// Convex hull
TEST(Geometry, ConvexHull)
{
    typedef int V;
    typedef Vec2D<V> T;
    const std::vector<T> p = {T(1,0),T(1,1),T(0,1),T(-1,1),T(-1,0),T(-1,-1),T(0,-1),T(1,-1)};
    const ConvexHull<V> ch(p);
    const std::vector<T> q = {T(-1,-1),T(1,-1),T(1,1),T(-1,1)};
    EXPECT_EQ(q, ch());
    EXPECT_EQ(4, ch.area());
}

// Line segment intersection
TEST(Geometry, LineSegmentIntersection)
{
    typedef Line2D<int> T;
    const T a(1,1,4,4);
    EXPECT_TRUE (intersect(a,T(1,2,2,1)));
    EXPECT_TRUE (intersect(a,T(3,2,2,3)));
    EXPECT_TRUE (intersect(a,T(3,4,4,3)));
    EXPECT_TRUE (intersect(a,T(0,0,1,1)));
    EXPECT_TRUE (intersect(a,T(2,2,1,1)));
    EXPECT_TRUE (intersect(a,T(2,2,3,3)));
    EXPECT_FALSE(intersect(a,T(5,5,6,6)));
    EXPECT_FALSE(intersect(a,T(0,0,-1,-1)));
    EXPECT_FALSE(intersect(a,T(0,0,1,-1)));
    EXPECT_FALSE(intersect(a,T(0,0,-1,1)));
    EXPECT_TRUE (intersect(a,T(0,2,1,1)));
    EXPECT_TRUE (intersect(a,T(2,0,1,1)));
    EXPECT_TRUE (intersect(a,T(0,2,2,0)));
    EXPECT_FALSE(intersect(a,T(0,1,1,0)));
    EXPECT_FALSE(intersect(a,T(0,2,-1,3)));
    EXPECT_FALSE(intersect(a,T(2,0,3,-1)));
    EXPECT_FALSE(intersect(a,T(1,2,3,4)));
    EXPECT_FALSE(intersect(a,T(2,1,4,3)));
    EXPECT_FALSE(intersect(a,T(1,4,2,3)));
    EXPECT_FALSE(intersect(a,T(3,2,4,1)));
}

// Triangle
TEST(Geometry, Triangle)
{
    EXPECT_EQ(0.5, Triangle2D<double>(1,1,0,0,1,0).area());
    EXPECT_EQ(  6, Triangle2D<double>(3,4,3,0,0,4).area());
    const auto EPS = 1e-9;
    const Triangle2D<double> t(0,0,2,0,2,2);
    auto p = t.circumcenter();
    EXPECT_NEAR(1.0, p.x, EPS);
    EXPECT_NEAR(1.0, p.y, EPS);
    EXPECT_NEAR(std::sqrt(2), t.R(), EPS);
    EXPECT_TRUE (t.inside(Vec2D<double>( 1.5, 0.5)));
    EXPECT_TRUE (t.inside(Vec2D<double>( 0.0, 0.0)));
    EXPECT_FALSE(t.inside(Vec2D<double>( 1.0, 1.5)));
    EXPECT_FALSE(t.inside(Vec2D<double>( 2.0, 3.0)));
    EXPECT_FALSE(t.inside(Vec2D<double>( 3.0, 1.0)));
    EXPECT_FALSE(t.inside(Vec2D<double>( 3.0,-1.0)));
    EXPECT_FALSE(t.inside(Vec2D<double>( 1.0,-0.5)));
    EXPECT_FALSE(t.inside(Vec2D<double>(-1.0,-0.5)));
}

// Probability of Complete Gacha
TEST(Probability, CompleteGacha)
{
    std::vector<int> a = { 1, 1 };
    EXPECT_EQ(3, comp_gacha(a));
    EXPECT_EQ(3, comp_gacha_avg(2));
    std::vector<int> b = { 1, 4, 5 };
    auto s145 = comp_gacha(b);
    EXPECT_LT(10.722222, s145);EXPECT_LT(s145, 10.722223);
    auto avg50 = comp_gacha_avg(50);
    EXPECT_LT(224.960266, avg50);EXPECT_LT(avg50, 224.960267);
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
