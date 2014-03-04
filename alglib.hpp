#include <vector>
#include <algorithm>
#include <functional>
#include <numeric>
#include <cmath>

// Greatest Common Divisor
template<class T>
T gcd(T a, T b)
{
    while(b)
        std::swap(a%=b,b);
    return a;
}
// Least Common Multiple
template<class T>
T lcm(T a, T b)
{
    return a/gcd(a,b)*b;
}

// Extended Greatest Common Divisor
template<class T>
T extgcd(T a, T b, T& x, T& y)
{
    T s=0, t=1;
    T u=1, v=0;
    while(a)
    {
        T w=b/a;
        std::swap(b-=w*a,a);
        std::swap(s-=w*t,t);
        std::swap(u-=w*v,v);
    }
    x=s;
    y=u;
    return b;
}

// Prime
template<class T>
bool prime(T n)
{
    for(T i=2; i*i<=n; ++i)
        if(n%i==0)
            return false;
    return true;
}

// Prime Factors
template<class T>
std::vector<T> primes(T n)
{
    std::vector<T> ps;
    for(T i=2; i*i<=n; ++i)
        while(n%i==0)
            n/=i,ps.push_back(i);
    if(n>1) ps.push_back(n);
    return ps;
}

// Euler's phi function
template<class T>
T euler_phi(T n)
{
    T p=n;
    for(T i=2; i*i<=n; ++i)
        if(n%i==0)
        {
            p=p/i*(i-1);
            do n/=i; while(n%i==0);
        }
    if(n>1) p=p/n*(n-1);
    return p;
}

// Integer Square Root
template<class T>
T isqrt(const T& i)
{
    if(i<=0) return 0;
    T x=1, y=i;
    while(x<=y) x<<=1, y>>=1;
    do std::swap(x, y=(x+i/x)>>1);
    while(x<y);
    return y;
}

// Power
template<class T, class U>
T POW(T a, U n)
{
    T t=1;
    while(n)
    {
        if(n&1) t*=a;
        a*=a;
        n>>=1;
    }
    return t;
}
// Combination
template<class T, class U=T>
U C(T n, T k)
{
    U r=1;
    k=std::max<T>(k,n-k);
    for(T i=n; i>k; --i)
        r*=i;
    k=n-k;
    for(T i=2; i<=k; ++i)
        r/=i;
    return r;
}

// Modulo Integer
template<class T1, class T2, bool B=(sizeof(T1)>=sizeof(T2))> struct max_type{ typedef T1 type; };
template<class T1, class T2> struct max_type<T1,T2,false>{ typedef T2 type; };
template<class T=int, T M=1000000007, class U=typename max_type<T,long long>::type>
class mint
{
    T v;
    struct Mul{ T v; Mul(T a, T b):v(U(a)*b%M){} };
    mint(const Mul& m):v(m.v){}
    inline static mint add(T a, T b){ return a+b; }
    inline static mint sub(T a, T b){ return add(a,M-b); }
    inline static mint mul(T a, T b){ return Mul(a,b); }
    inline static mint div(T a, T b){ return mul(a,inv(b)); }
    // gcd(a,M) must be 1
    static T inv(T a)
    {
        T b=M;
        T u=0, v=1;
        while(a)
        {
            T t=b/a;
            std::swap(b-=t*a,a);
            std::swap(u-=t*v,v);
        }
        return (u+M)%M;
    }
public:
    typedef T value_type;
    T operator()() const { return v; }
    mint(const mint& m):v(m.v){}
    mint(T i=T()):v(i%M){}
    mint& operator+=(const mint& m){ return *this = add(v,m.v); }
    mint& operator-=(const mint& m){ return *this = sub(v,m.v); }
    mint& operator*=(const mint& m){ return *this = mul(v,m.v); }
    mint& operator/=(const mint& m){ return *this = div(v,m.v); }
    friend mint operator+(const mint& l, const mint& r){ return add(l.v,r.v); }
    friend mint operator-(const mint& l, const mint& r){ return sub(l.v,r.v); }
    friend mint operator*(const mint& l, const mint& r){ return mul(l.v,r.v); }
    friend mint operator/(const mint& l, const mint& r){ return div(l.v,r.v); }
    mint pow(T n) const { return POW(*this,n); }
    static mint pow(T a, T n){ return POW<mint,T>(a,n); }
    static mint c(T n, T k){ return C<T,mint>(n,k); }
};

// Longest Common Subsequence
template<class T, class S=int>
S lcs(const T& a, const T& b)
{
    const int N = b.size();
    std::vector<S> dp(N+1,0);
    std::vector<S> dpn(N+1,0);
    for(auto v: a)
    {
        for(int i=0; i<N; ++i)
            dpn[i+1]=std::max(dpn[i],dp[i+(v==b[i]?0:1)]+(v==b[i]?1:0));
        dp.swap(dpn);
    }
    return dp[N];
}

// Warshall-Floyd
template<class T> struct minimum{ const T& operator()(const T& a, const T& b){ return std::min<T>(a,b); } };
template<template<class> class F=std::plus, template<class> class G=minimum, class T, class U=typename T::value_type::value_type>
void wf(T& g)
{
    const int N = g.size();
    for(int k=0; k<N; ++k)
        for(int i=0; i<N; ++i)
            for(int j=0; j<N; ++j)
                g[i][j] = G<U>()(g[i][j],F<U>()(g[i][k],g[k][j]));
}

// Area of triangle
double area_of_triangle(double x1, double y1, double x2, double y2, double x3, double y3)
{
    return std::abs(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2))/2;
}

// Probability of Complete Gacha
template<class T=double, class U>
T comp_gacha(U ns)
{
    const int N = ns.size();
    const T s = std::accumulate(ns.begin(),ns.end(),T());
    std::vector<T> p(N);
    std::transform(ns.begin(),ns.end(),p.begin(),[&s](const T& x){ return x/s; });
    std::vector<T> dp(1<<N);
    for(int b=1; b<1<<N; ++b)
    {
        T t=0, r=0;
        for(int i=0; i<N; ++i)
            if(b&1<<i)
                r+=p[i], t+=dp[b^1<<i]*p[i];
        dp[b]=(1+t)/r;
    }
    return dp[(1<<N)-1];
}

template<class T=double>
T comp_gacha_avg(int n)
{
    T t=0;
    for(int i=n; i>0; --i)
        t += static_cast<T>(n)/i;
    return t;
}
