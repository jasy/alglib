#include <vector>
#include <algorithm>

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
