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
    // Prime
    assert( prime(2));
    assert( prime(3));
    assert(!prime(4));
    assert( prime(5));
    assert(!prime(6));
    assert( prime(7));
    assert(!prime(8));
    assert(!prime(9));
    assert( prime(9999991));
    assert(!prime(99999991));
    assert( prime(129402307));
    assert( prime(282437925089LL));
    // Prime Factors
    assert(primes(103)==std::vector<int>({103}));
    assert(primes(2*3*3*7*13*13*13)==std::vector<int>({2,3,3,7,13,13,13}));
    // Euler's phi function
    assert(euler_phi(1)==1);
    assert(euler_phi(2)==1);
    assert(euler_phi(3)==2);
    assert(euler_phi(4)==2);
    assert(euler_phi(5)==4);
    assert(euler_phi(6)==2);
    assert(euler_phi(7)==6);
    assert(euler_phi(8)==4);
    assert(euler_phi(9)==6);
    assert(euler_phi(10)==4);
    assert(euler_phi(11)==10);
    assert(euler_phi(12)==4);
    assert(euler_phi(13)==12);
    assert(euler_phi(14)==6);
    assert(euler_phi(15)==8);
    assert(euler_phi(16)==8);
    assert(euler_phi(5186)==2592);
    assert(euler_phi(5187)==2592);
    assert(euler_phi(5188)==2592);
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
    // Square Matrix
    {
        enum{ N=3 };
        SquareMatrix<int> m(N);
        for(int i=0; i<N; ++i)
            m(i,i) = 1;
        m*=m*m;
        for(int r=0; r<N; ++r) for(int c=0; c<N; ++c)
            { assert(m(r,c)==(r==c?1:0)); }
        {
            auto t = m.pow(1000000000);
            for(int r=0; r<N; ++r) for(int c=0; c<N; ++c)
                { assert(t(r,c)==(r==c?1:0)); }
        }
        for(int r=0; r<N; ++r) for(int c=0; c<N; ++c)
            m(r,c) = r+1;
        {
            auto t = m*m;
            for(int r=0; r<N; ++r) for(int c=0; c<N; ++c)
                { assert(t(r,c)==N*(N+1)/2*(r+1)); }
        }
        {
            auto t = m.pow(2);
            for(int r=0; r<N; ++r) for(int c=0; c<N; ++c)
                { assert(t(r,c)==N*(N+1)/2*(r+1)); }
        }
        std::vector<int> v(N);
        for(int i=0; i<N; ++i) v[i]=i+1;
        {
            auto t = m*v;
            for(int i=0; i<N; ++i)
                { assert(t[i]==N*(N+1)/2*(i+1)); }
        }
    }
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
        assert(2148==(m.pow(999999999-N)*A)[N-1]);
    }
    // Fenwick Tree (Binary Indexed Tree)
    {
        enum{ N=10 };
        fenwick<int> t(N);
        t.add(0,1);
        t.add(N/2,1);
        assert(0==t.sum(0));
        assert(1==t.sum(1));
        assert(1==t.sum(N/2));
        assert(2==t.sum(N/2+1));
        assert(2==t.sum(N));
        assert(0==t.lower_bound(0));
        assert(0==t.lower_bound(1));
        assert(N/2==t.lower_bound(2));
        assert(N==t.lower_bound(3));
        t.add(N/2,-1);
        t.add(0,-1);
        for(int i=0; i<N; ++i) t.add(i,i+1);
        assert(N*(N+1)/2==t.sum(N));
        for(int i=0; i<N; ++i){ assert(i*(i+1)/2==t.sum(i)); }
        for(int i=0; i<N; ++i) t.add(i,N-(i+1));
        assert(N*N==t.sum(0,N));
        for(int i=0; i<N; ++i){ assert(N*(N-i)==t.sum(i,N)); }
    }
    {
        enum{ N=10 };
        fenwick_range<int> t(N);
        for(int i=0; i<N; ++i) t.add(i,N,1);
        assert(N*(N+1)/2==t.sum(N));
        for(int i=0; i<N; ++i){ assert(i*(i+1)/2==t.sum(i)); }
        for(int i=0; i<N; ++i) t.add(0,i,1);
        assert(N*N==t.sum(0,N));
        for(int i=0; i<N; ++i){ assert(N*(N-i)==t.sum(i,N)); }
    }
    {
        enum{ N=3, M=5 };
        fenwick2D<int> t(N,M);
        for(int i=0; i<N; ++i) for(int j=0; j<M; ++j) t.add(i,j,1);
        assert(N*M==t.sum(N,M));
        for(int i=0; i<N; ++i) for(int j=0; j<M; ++j){ assert(i*j==t.sum(i,j)); }
        for(int i1=0; i1<N; ++i1) for(int j1=0; j1<M; ++j1)
            for(int i2=i1; i2<N; ++i2) for(int j2=j1; j2<M; ++j2)
                { assert((i2-i1)*(j2-j1)==t.sum(i1,j1,i2,j2)); }
    }
    // Longest Common Subsequence
    assert(3==lcs<std::string>("13579","395678"));
    {
        std::vector<int> a = {  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27 };
        std::vector<int> b = {  2,  3, 30,  4,  5,  6, 31,  7, 32,  9, 10, 11, 12, 33, 14, 15, 16, 34, 17, 19, 20, 22, 35, 24, 36, 37, 25, 26, 27, 38, 39 };
        assert(21==lcs(a,b));
    }
    // Longest Increasing Subsequence
    assert(4==lis(std::string("135246")));
    assert(1==lis(std::string("54321")));
    assert(7==lis(std::string("1234567")));
    assert(6==lis(std::vector<int>({1,3,4,5,6,7,2})));
    // Warshall-Floyd
    {
        const int N = 10;
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
    // Shortest Path Faster Algorithm
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
        assert(N==cost.size());
        for(int i=0; i<N; ++i){ assert(N-1-i==cost[i]); }
    }
    // Dijkstra's algorithm
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
        assert(N==cost.size());
        for(int i=0; i<N; ++i){ assert(N-1-i==cost[i]); }
    }
    // Lowest Common Ancestor
    {
        std::vector<std::pair<int,int>> edges1 = {{2,0},{1,0},{1,3},{1,4},{2,5},{2,6}};
        const LCA lca1(edges1);
        assert(2==lca1.distance(3,4));
        assert(2==lca1.distance(0,5));
        assert(4==lca1.distance(4,5));
        assert(4==lca1.distance(3,6));
        assert(3==lca1.distance(4,2));
        std::vector<std::pair<int,int>> edges2 = {{0,1},{0,2},{0,3},{3,4}};
        const LCA lca2(edges2);
        assert(2==lca2.distance(1,2));
        assert(2==lca2.distance(1,3));
        assert(3==lca2.distance(1,4));
    }
    // Union Find
    {
        {
            UnionFind<int> uf(6);
            assert(!uf.same(2,5));
            std::vector<std::pair<int,int>> edges = {{0,4},{0,5},{1,3},{1,4},{2,3}};
            for(auto& e: edges)
                uf.unite(e.first,e.second);
            for(auto& e: edges)
                assert(uf.same(e.first,e.second));
            assert( uf.same(2,5));
        }
        {
            UnionFind<int> uf(6);
            assert(!uf.same(2,5));
            assert(!uf.same(1,2));
            assert(!uf.same(4,5));
            std::vector<std::pair<int,int>> edges = {{0,4},{0,5},{1,3},{2,3}};
            for(auto& e: edges)
                uf.unite(e.first,e.second);
            for(auto& e: edges)
                assert(uf.same(e.first,e.second));
            assert(!uf.same(2,5));
            assert( uf.same(1,2));
            assert( uf.same(4,5));
        }
    }
    // Kruskal's algorithm
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
        assert(350==t);
    }
    // Dinic's algorithm
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
            assert(1==dinic(edges,N,0,N-1));
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
            assert(1==dinic(edges,N,0,N-1));
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
            assert(2==dinic(edges,N,0,N-1));
        }
    }
    // Bipartite Matching
    {
        std::vector<std::pair<int,int>> edges1 = {{0,3},{0,4},{0,5},{1,3},{1,4},{2,3}};
        assert(3==bm(edges1,6));
        std::multimap<int,int> edges2 = {{0,4},{0,6},{0,7},{1,4},{1,5},{1,8},{2,5},{2,8},{3,6}};
        assert(4==bm(edges2,9));
        typedef std::tuple<int,int> T3;
        std::array<T3,1> edges3 = {{T3{0,1}}};
        assert(1==bm(edges3,2));
    }
    // Line segment intersection
    {
        typedef Line2D<int> T;
        const T a(1,1,4,4);
        assert( intersect(a,T(1,2,2,1)));
        assert( intersect(a,T(3,2,2,3)));
        assert( intersect(a,T(3,4,4,3)));
        assert( intersect(a,T(0,0,1,1)));
        assert( intersect(a,T(2,2,1,1)));
        assert( intersect(a,T(2,2,3,3)));
        assert(!intersect(a,T(5,5,6,6)));
        assert(!intersect(a,T(0,0,-1,-1)));
        assert(!intersect(a,T(0,0,1,-1)));
        assert(!intersect(a,T(0,0,-1,1)));
        assert( intersect(a,T(0,2,1,1)));
        assert( intersect(a,T(2,0,1,1)));
        assert( intersect(a,T(0,2,2,0)));
        assert(!intersect(a,T(0,1,1,0)));
        assert(!intersect(a,T(0,2,-1,3)));
        assert(!intersect(a,T(2,0,3,-1)));
        assert(!intersect(a,T(1,2,3,4)));
        assert(!intersect(a,T(2,1,4,3)));
        assert(!intersect(a,T(1,4,2,3)));
        assert(!intersect(a,T(3,2,4,1)));
    }
    // Triangle
    {
        assert(0.5==Triangle2D<double>(1,1,0,0,1,0).area());
        assert(  6==Triangle2D<double>(3,4,3,0,0,4).area());
        const auto EPS = 1e-9;
        const Triangle2D<double> t(0,0,2,0,2,2);
        auto p = t.circumcenter();
        assert(1.0-EPS<=p.x && p.x<=1.0+EPS);
        assert(1.0-EPS<=p.y && p.y<=1.0+EPS);
        auto R = t.R();
        auto sqrt2 = std::sqrt(2);
        assert(sqrt2-EPS<=R && R<=sqrt2+EPS);
        assert( t.inside(Vec2D<double>( 1.5, 0.5)));
        assert( t.inside(Vec2D<double>( 0.0, 0.0)));
        assert(!t.inside(Vec2D<double>( 1.0, 1.5)));
        assert(!t.inside(Vec2D<double>( 2.0, 3.0)));
        assert(!t.inside(Vec2D<double>( 3.0, 1.0)));
        assert(!t.inside(Vec2D<double>( 3.0,-1.0)));
        assert(!t.inside(Vec2D<double>( 1.0,-0.5)));
        assert(!t.inside(Vec2D<double>(-1.0,-0.5)));
    }
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
