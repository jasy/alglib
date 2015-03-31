#include <vector>
#include <queue>
#include <algorithm>
#include <functional>
#include <numeric>
#include <limits>
#include <cmath>
#include <type_traits>

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

// Divisors
template<class T>
std::vector<T> divisors(T n)
{
    std::vector<T> ds,lds;
    for(T i=1; i*i<=n; ++i)
        if(n%i==0)
            ds.push_back(i),
            lds.push_back(n/i);
    const int m = ds[ds.size()-1]==lds[lds.size()-1]?1:0;
    ds.reserve(ds.size()*2-m);
    ds.insert(ds.end(),lds.rbegin()+m,lds.rend());
    return ds;
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
template<class T1, class T2> struct max_type{ typedef typename std::conditional<sizeof(T1)>=sizeof(T2),T1,T2>::type type; };
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

// Square Matrix
template<class T, template<class> class Plus=std::plus, template<class> class Mul=std::multiplies, T PlusUnit=0, T MulUnit=1>
class SquareMatrix
{
    std::vector<T> a;
    static void update(T& a, T b, T c){ a = Plus<T>()(a, Mul<T>()(b,c)); }
    SquareMatrix& operator=(SquareMatrix&& m){ this->a = std::move(m.a); return *this; }
public:
    const size_t N;
    explicit SquareMatrix(size_t N):a(N*N,PlusUnit),N(N){}
    SquareMatrix(const SquareMatrix& m):a(m.a),N(m.N){}
    SquareMatrix(SquareMatrix&& m):a(std::move(m.a)),N(m.N){}
    T& operator()(size_t r, size_t c){ return a[r*N+c]; }
    T operator()(size_t r, size_t c)const{ return a[r*N+c]; }
    std::vector<T> operator*(const std::vector<T>& v)const
    {
        const auto N=this->N;
        assert(N==v.size());
        auto& m = *this;
        std::vector<T> v2(N);
        for(size_t r=0; r<N; ++r)
            for(size_t c=0; c<N; ++c)
                update(v2[r],m(r,c),v[c]);
        return v2;
    }
    SquareMatrix& operator*=(const SquareMatrix& m){ return *this = std::move(*this*m); }
    friend SquareMatrix operator*(const SquareMatrix& lhs, const SquareMatrix& rhs)
    {
        const auto N=lhs.N;
        assert(N==rhs.N);
        SquareMatrix m(N);
        for(size_t r=0; r<N; ++r)
            for(size_t c=0; c<N; ++c)
                for(size_t i=0; i<N; ++i)
                    update(m(r,c),lhs(r,i),rhs(i,c));
        return m;
    }
    SquareMatrix pow(unsigned int n)const
    {
        SquareMatrix a(*this);
        SquareMatrix t(N); for(size_t i=0; i<N; ++i) t(i,i)=MulUnit;
        while(n)
        {
            if(n&1) t*=a;
            a*=a;
            n>>=1;
        }
        return t;
    }
};

// Fenwick Tree (Binary Indexed Tree)
#ifdef _MSC_VER
#pragma warning(disable:4146)
#endif
template<class T>
class fenwick
{
    std::vector<T> a;
    static size_t msb(size_t v)
    {
        for(size_t i=1; i<sizeof(size_t)*8; i<<=1) v |= v>>i;
        return (v>>1)+1;
    }
public:
    explicit fenwick(size_t N):a(N+1,0){}
    void add(size_t i, T v)
    {
        const size_t N = a.size()-1;
        for(size_t x=i+1; x<=N; x+=x&-x) a[x]+=v;
    }
    T sum(size_t i)
    {
        T t=0;
        for(size_t x=i; x>0; x-=x&-x) t+=a[x];
        return t;
    }
    T sum(size_t i, size_t j){ return sum(j)-sum(i); }
    size_t lower_bound(T v)
    {
        if(v<=0) return 0;
        const size_t N = a.size()-1;
        size_t i=0;
        for(size_t k=msb(N); k>0; k>>=1)
            if (i+k<=N && v>a[i+k])
            {
                v -= a[i+k];
                i += k;
            }
        return i;
    }
};
template<class T>
class fenwick_range
{
    fenwick<T> p,q;
public:
    explicit fenwick_range(size_t N):p(N+1),q(N+1){}
    void add(size_t i, size_t j, T v)
    {
        p.add(i,-v*T(i));
        p.add(j, v*T(j));
        q.add(i, v);
        q.add(j,-v);
    }
    T sum(size_t i){ return p.sum(i+1)+q.sum(i+1)*T(i); }
    T sum(size_t i, size_t j){ return sum(j)-sum(i); }
};
template<class T>
class fenwick2D
{
    size_t N_,M_;
    std::vector<T> a;
public:
    explicit fenwick2D(size_t N, size_t M):N_(N),M_(M),a((N+1)*(M+1),0){}
    void add(size_t i, size_t j, T v)
    {
        const size_t N=N_, M=M_;
        for(size_t x=i+1; x<=N; x+=x&-x)
            for(size_t y=j+1; y<=M; y+=y&-y)
                a[x*(M+1)+y]+=v;
    }
    T sum(size_t i, size_t j)
    {
        T t=0;
        const size_t M=M_;
        for(size_t x=i; x>0; x-=x&-x)
            for(size_t y=j; y>0; y-=y&-y)
                t+=a[x*(M+1)+y];
        return t;
    }
    T sum(size_t i1, size_t j1, size_t i2, size_t j2)
    {
        return sum(i2,j2)-sum(i2,j1)-sum(i1,j2)+sum(i1,j1);
    }
};
#ifdef _MSC_VER
#pragma warning(default:4146)
#endif

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

// Longest Increasing Subsequence
template<class T, class S=typename T::value_type>
int lis(const T& a)
{
    const S INF = std::numeric_limits<S>::max();
    std::vector<S> dp(a.size(),INF);
    for(const auto v: a)
        *std::lower_bound(dp.begin(),dp.end(),v) = v;
    return std::lower_bound(dp.begin(),dp.end(),INF)-dp.begin();
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

// Shortest Path Faster Algorithm
template<class Cost=int, class Vertex=int, class E=std::pair<Cost,Vertex>>
struct SPFA
{
    typedef E Edge;
    static std::vector<Cost> spfa(const std::vector<std::vector<Edge>>& edges, Vertex s)
    {
        const Cost INF = std::numeric_limits<Cost>::max();
        const Vertex N = edges.size();
        std::vector<Cost> d(N,INF);
        std::queue<Vertex> q;
        std::vector<bool> r(N,false);
        d[s]=0; r[s]=true; q.push(s);
        while(!q.empty())
        {
            Vertex u = q.front();q.pop();
            r[u]=false;
            for(const Edge&e:edges[u])
            {
                const Cost& c = std::get<0>(e);
                const Vertex& v = std::get<1>(e);
                if(d[v] > d[u]+c)
                {
                    d[v] = d[u]+c;
                    if(!r[v])
                    {
                        r[v]=true;
                        q.push(v);
                    }
                }
            }
        }
        return d;
    }
};

// Dijkstra's algorithm
template<class Cost=int, class Vertex=int, class E=std::pair<Cost,Vertex>>
struct Dijkstra
{
    typedef E Edge;
    static std::vector<Cost> dijkstra(const std::vector<std::vector<Edge>>& edges, Vertex s)
    {
        const Cost INF = std::numeric_limits<Cost>::max();
        const Vertex N = edges.size();
        std::vector<Cost> d(N,INF);
        std::priority_queue<Edge,std::vector<Edge>,std::greater<Edge>> q;
        q.emplace(0,s);
        while(!q.empty())
        {
            const Edge& t = q.top();
            const Cost b = std::get<0>(t);
            const Vertex u = std::get<1>(t);
            q.pop();
            if(d[u]!=INF) continue;
            d[u]=b;
            for(const Edge&e:edges[u])
            {
                const Cost& c = std::get<0>(e);
                const Vertex& v = std::get<1>(e);
                if(d[v]==INF)
                    q.emplace(b+c,v);
            }
        }
        return d;
    }
};

// Lowest Common Ancestor
class LCA
{
    typedef int V;
    const V N;
    const int l2;
    std::vector<std::vector<V>> g;
    std::vector<int> r;
    std::vector<std::vector<V>> p;
    void make(const V s, const V t, const int d)
    {
        r[t]=d;
        for(auto v: g[t]) if(v!=s)
        {
            p[v][0]=t;
            make(t,v,d+1);
        }
    }
    static int lg(const V n)
    {
        int i=0;
        while(n>(V(1)<<i)+1) ++i;
        return i+1;
    }
    LCA():N(0),l2(0){}
public:
    LCA(const std::vector<std::pair<V,V>>& edges)
        :N(edges.size()+1),l2(lg(N)),
        g(N),r(N),p(N,std::vector<int>(l2))
    {
        for(auto& v: edges)
        {
            V x,y; std::tie(x,y)=v;
            g[x].emplace_back(y);
            g[y].emplace_back(x);
        }
        make(N,0,0);
        for(int i=0; i<l2-1; ++i)
            for(V j=0; j<N; ++j)
                p[j][i+1] = p[p[j][i]][i];
    }
    V lca(V a, V b) const
    {
        if(r[a]<r[b]) std::swap(a,b);
        for(int i=0; i<l2; ++i)
            if((r[a]-r[b]) & (1<<i)) a=p[a][i];
        if(a==b) return a;
        for(int i=l2-1; i>=0; --i)
            if(p[a][i]!=p[b][i]) a=p[a][i],b=p[b][i];
        return p[a][0];
    }
    int distance(const V a, const V b) const
    {
        return r[a]+r[b]-2*r[lca(a,b)];
    }
};

// Union Find
template<class T=int>
class UnionFind
{
    std::vector<T> p;
    std::vector<T> r;
    T find(const T x)
    {
        if(x==p[x]) return x;
        return p[x]=find(p[x]);
    }
    void unite_(T x, T y)
    {
        if(x==y) return;
        if(r[x]>r[y]) std::swap(x,y);
        else if(r[x]==r[y]) ++r[y];
        p[x]=y;
    }
public:
    explicit UnionFind(const T n):p(n),r(n,0)
    {
        for(T i=0; i<n; ++i)
            p[i]=i;
    }
    bool same(const T x, const T y)
    {
        // assert(0<=x && x<p.size())
        // assert(0<=y && y<p.size())
        return find(x)==find(y);
    }
    void unite(const T x, const T y)
    {
        // assert(0<=x && x<p.size())
        // assert(0<=y && y<p.size())
        unite_(find(x),find(y));
    }
};

// Kruskal's algorithm
template<class T, class U>
T& kruskal(T& edges, U max_v)
{
    typedef typename T::value_type Edge;
    std::sort(edges.begin(),edges.end());
    UnionFind<U> uf(max_v);
    edges.erase(
        std::remove_if(edges.begin(),edges.end(),
            [&uf](const Edge& e){
                const auto& a = std::get<1>(e);
                const auto& b = std::get<2>(e);
                if(uf.same(a,b)) return true;
                uf.unite(a,b);
                return false;
            }),
        edges.end());
    return edges;
}

// Dinic's algorithm
template<class T, class U=int, class F=int>
F dinic(const T& edges, U max_v, U s, U t)
{
    typedef typename T::value_type Edge;
    std::vector<std::vector<Edge>> g(max_v);
    for(auto& e: edges){
        const auto& c = std::get<0>(e);
        const auto& s = std::get<1>(e);
        const auto& t = std::get<2>(e);
        g[s].emplace_back(c, t, g[t].size());
        g[t].emplace_back(0, s, g[s].size()-1);
    }
    const F INF = std::numeric_limits<F>::max()/2;
    F flow=0;
    for(;;){
        std::vector<int> dist(max_v,-1); dist[s]=0;
        {
            std::queue<U> q; q.push(s);
            while(!q.empty()){
                auto v = q.front(); q.pop();
                for(auto& e: g[v]){
                    const auto& c = std::get<0>(e);
                    const auto& t = std::get<1>(e);
                    if(c>0 && dist[t]<0)
                        dist[t] = dist[v]+1,
                        q.push(t);
                }
            }
        }
        if(dist[t]<0) break;
        std::vector<int> iter(max_v,0);
        class Func{
        public:
            static F dfs(U v, U u, F f, std::vector<std::vector<Edge>>& g, const std::vector<int>& dist, std::vector<int>& iter){
                if(v==u) return f;
                const int N = g[v].size();
                for(auto& i=iter[v]; i<N; ++i){
                    auto& e = g[v][i];
                    auto& c = std::get<0>(e);
                    const auto& t = std::get<1>(e);
                    const auto& r = std::get<2>(e);
                    if(c>0 && dist[v] < dist[t]){
                        F d = dfs(t, u, std::min(f,c), g, dist, iter);
                        if(d!=0){
                            c-=d;
                            std::get<0>(g[t][r])+=d;
                            return d;
                        }
                    }
                }
                return 0;
            }
        };
        F f;
        while((f = Func::dfs(s, t, INF, g, dist, iter))>0) flow += f;
    }
    return flow;
}

// Bipartite Matching
template<class T, class U>
int bm(const T& edges, U max_v)
{
    std::vector<std::vector<U>> g(max_v,std::vector<U>());
    std::vector<U> match(max_v,max_v);
    std::vector<bool> used;
    std::function<bool(U)> dfs;
    dfs = [&](U v){
        used[v]=true;
        for(auto u: g[v]){
            auto w = match[u];
            if(w==max_v || (!used[w] && dfs(w))){
                match[v]=u;
                match[u]=v;
                return true;
            }
        }
        return false;
    };
    for(auto& p: edges)
    {
        U a,b; std::tie(a,b)=p;
        g[a].push_back(b),
        g[b].push_back(a);
    }
    int c=0;
    for(U v=0; v<max_v; ++v)
        if(match[v]==max_v)
        {
            used.assign(max_v,false);
            if(dfs(v))
                ++c;
        }
    return c;
}

// Vector 2D
template<class T>
struct Vec2D
{
    T x,y;
    explicit Vec2D(T x_=T(), T y_=T()):x(x_),y(y_){}
    Vec2D(const std::pair<T,T>& p):x(p.first),y(p.second){}
    Vec2D& operator+=(const Vec2D& v){ return *this = *this+v; }
    Vec2D& operator-=(const Vec2D& v){ return *this = *this-v; }
    friend Vec2D operator+(const Vec2D& l, const Vec2D& r){ return Vec2D(l.x+r.x,l.y+r.y); }
    friend Vec2D operator-(const Vec2D& l, const Vec2D& r){ return Vec2D(l.x-r.x,l.y-r.y); }
    friend T abs(const Vec2D& v){ return std::hypot(v.x,v.y); }
    friend T norm(const Vec2D& v){ return v.x*v.x+v.y*v.y; }
    friend T dot(const Vec2D& l, const Vec2D& r){ return l.x*r.x+l.y*r.y; }
    friend T cross(const Vec2D& l, const Vec2D& r){ return l.x*r.y-l.y*r.x; }
    friend int dir(const Vec2D& l, const Vec2D& r)
    {
        // assert(norm(l)!=0);
        const auto c = cross(l,r);
        if(c>0) return 1;
        if(c<0) return -1;
        if(dot(l,r)<0) return -2;
        if(norm(l)<norm(r)) return 2;
        return 0;
    }
};

// Line segment 2D
template<class T>
class Line2D
{
    Vec2D<T> a,b;
    static bool judge(const Line2D& l, const Line2D& r)
    {
        auto v = l.b-l.a;
        return dir(v,r.a-l.a)*dir(v,r.b-l.a)<=0;
    }
public:
    explicit Line2D(T ax=T(), T ay=T(), T bx=T(), T by=T()):a(ax,ay),b(bx,by){}
    Line2D(const Vec2D<T>& a_, const Vec2D<T>& b_):a(a_),b(b_){}
    friend bool intersect(const Line2D& l, const Line2D& r)
    {
        return judge(l,r) && judge(r,l);
    }
};

// Triangle 2D
template<class T=double>
class Triangle2D
{
    Vec2D<T> a,b,c;
public:
    explicit Triangle2D(T ax=T(), T ay=T(), T bx=T(), T by=T(), T cx=T(), T cy=T()):a(ax,ay),b(bx,by),c(cx,cy){}
    Triangle2D(const Vec2D<T>& a_, const Vec2D<T>& b_, const Vec2D<T>& c_):a(a_),b(b_),c(c_){}
    T area() const
    {
        return std::abs(a.x*(b.y-c.y)+b.x*(c.y-a.y)+c.x*(a.y-b.y))/2;
    }
    T dist_a() const { return abs(b-c); }
    T dist_b() const { return abs(c-a); }
    T dist_c() const { return abs(a-b); }
    Vec2D<T> circumcenter() const
    {
        const Vec2D<T> vb(2*(b.x-a.x),2*(c.x-a.x)),
                       vc(2*(b.y-a.y),2*(c.y-a.y)),
                       vd(norm(a)-norm(b),norm(a)-norm(c));
        return Vec2D<T>(cross(vc,vd)/cross(vb,vc),
                        cross(vd,vb)/cross(vb,vc));
    }
    T R() const { return dist_a()*dist_b()*dist_c()/(4*area()); }
    bool inside(const Vec2D<T>& p, bool line=true) const
    {
        const auto c1 = cross(b-a,p-a);
        const auto c2 = cross(c-b,p-b);
        const auto c3 = cross(a-c,p-c);
        if(line) return (c1>=0&&c2>=0&&c3>=0)||(c1<=0&&c2<=0&&c3<=0);
        else     return (c1> 0&&c2> 0&&c3> 0)||(c1< 0&&c2< 0&&c3< 0);
    }
};

// Probability of Complete Gacha
template<class T=double, class U=std::vector<T>>
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
