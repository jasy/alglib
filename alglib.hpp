#include <vector>
#include <algorithm>

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
    k=std::max(k,n-k);
    for(T i=n; i>k; --i)
        r*=i;
    k=n-k;
    for(T i=2; i<=k; ++i)
        r/=i;
    return r;
}

// Modulo Integer
static constexpr int MOD = 1000000007;
template<class T1, class T2, bool B=(sizeof(T1)>=sizeof(T2))> struct max_type{ typedef T1 type; };
template<class T1, class T2> struct max_type<T1,T2,false>{ typedef T2 type; };
template<class T, T M=MOD, class U=typename max_type<T,long long>::type>
class mint
{
    T v;
    inline static T add(T a, T b){ return a+b; }
    inline static T sub(T a, T b){ return a-b+M; }
    inline static U mul(T a, T b){ return U(a)*b; }
    inline static U div(T a, T b){ return U(a)*inv(b); }
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
    T get() const { return v; }
    mint():v(0){}
    mint(const mint& m):v(m.v){}
    template<class S> mint(S i):v(i%M){}
    mint& operator+=(const mint& m){ return *this = add(v,m.v); }
    mint& operator-=(const mint& m){ return *this = sub(v,m.v); }
    mint& operator*=(const mint& m){ return *this = mul(v,m.v); }
    mint& operator/=(const mint& m){ return *this = div(v,m.v); }
    mint operator+(const mint& m) const { return add(v,m.v); }
    mint operator-(const mint& m) const { return sub(v,m.v); }
    mint operator*(const mint& m) const { return mul(v,m.v); }
    mint operator/(const mint& m) const { return div(v,m.v); }
    mint pow(T n) const { return POW(*this,n); }
    static mint pow(T a, T n){ return POW<mint,T>(a,n); }
    static mint c(T n, T k){ return C<T,mint>(n,k); }
};
template<class S, class T, T M, class U=T> mint<T,M,U> operator+(S l, mint<T,M,U> r){ return mint<T,M,U>(l)+=r; }
template<class S, class T, T M, class U=T> mint<T,M,U> operator-(S l, mint<T,M,U> r){ return mint<T,M,U>(l)-=r; }
template<class S, class T, T M, class U=T> mint<T,M,U> operator*(S l, mint<T,M,U> r){ return mint<T,M,U>(l)*=r; }
template<class S, class T, T M, class U=T> mint<T,M,U> operator/(S l, mint<T,M,U> r){ return mint<T,M,U>(l)/=r; }

// Longest Common Subsequence
template<class T, class S=int>
S lcs(const T& a, const T& b)
{
    const auto N = b.size();
    std::vector<S> dp(N+1,0);
    std::vector<S> dpn(N+1,0);
    for(auto v: a)
    {
        for(S i=0; i<N; ++i)
            dpn[i+1]=std::max(dpn[i],dp[i+(v==b[i]?0:1)]+(v==b[i]?1:0));
        dp.swap(dpn);
    }
    return dp[N];
}
