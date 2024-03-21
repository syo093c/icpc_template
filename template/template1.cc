#ifndef LOCAL
#pragma GCC optimize ("Ofast")
#pragma GCC optimize ("unroll-loops")
#endif

#include <bits/stdc++.h>
using namespace std;

using ll=long long;
#define rng(i,a,b) for(int i=int(a);i<int(b);i++)
#define rep(i,b) rng(i,0,b)
#define pb push_back
#define eb emplace_back
#define all(x) x.begin(),x.end()
#define si(x) int(x.size())
#define tct template<class t>
tct using vc=vector<t>;
#ifdef LOCAL
#define dmp(x) cerr<<__LINE__<<" "<<#x<<" "<<x<<endl
#else
#define dmp(x) void(0)
#endif

template<class t,class u> bool chmax(t&a, u b){ if(a<b){a=b; return true;} else return false;}
template<class t,class u> bool chmin(t&a, u b){ if(b<a){a=b; return true;} else return false;}

template<class t,class u>
ostream& operator<<(ostream& os, const pair<t,u>& p){
	return os<<"{"<<p.first<<","<<p.second<<"}";
}

template<class t> ostream& operator<<(ostream& os,const vector<t>& v){
	os<<"{";
	for(auto e: v)os<<e<<",";
	return os<<"}";
}

using P=pair<int,int>;
#define a first
#define b second
#define mp make_pair
using uint = unsigned int;
const uint mod=998244353;

struct mint{
	uint v;
	mint(ll vv=0){ s(vv%mod+mod); }
	mint& s(uint vv){
		v=vv<mod?vv:vv-mod;
		return *this;
	}
	mint operator-()const{ return mint()-*this;}
	mint&operator+=(mint r){ return s(v+r.v);}
	mint&operator-=(mint r){ return s(v+mod-r.v); }
	mint&operator*=(mint r){ v=(unsigned ll)v*r.v%mod; return *this;}
	mint&operator/=(mint r){ return *this*=r.inv();}
	mint operator+(mint r)const{return mint(*this)+=r;}
	mint operator-(mint r)const{return mint(*this)-=r;}
	mint operator*(mint r)const{return mint(*this)*=r;}
	mint operator/(mint r)const{return mint(*this)/=r;}
	mint pow(ll n)const{
		if(n<0)return inv().pow(-n);
		mint res(1),x(*this);
		while(n){
			if(n&1)res*=x;
			x*=x;
			n>>=1;
		}
		return res;
	}
	mint inv()const{return pow(mod-2);}
};
using ull = unsigned ll;
mt19937_64 mt(0);

int n,m;
pair<int,int>za[100005];
vc<P>edge[100005];;
vc<P>de[100005];
int ran[100005];
ull hsh[100005];
int deg[100005];
map<ull, pair<ll,int>>M;
map<P,int>SE;
void solve(){
	cin>>n>>m;
	rng(i,1,n+1) cin>>za[i].a>>za[i].b;
	rng(i,1,n+1){
		hsh[i] = mt();
	}
	ll ans = 0;
	rep(i,m){
		int u,v,w;cin>>u>>v>>w;
		SE[mp(u,v)] = w;
		SE[mp(v,u)] = w;
		edge[u].eb(v, w);
		edge[v].eb(u, w);
		deg[u] ++;
		deg[v] --;
		chmax(ans, w);
	}
	set<P>SS;
	rng(i,1,n+1) SS.insert(mp(deg[i], i));
	rng(r,1,n+1){
		auto [d,v] = *SS.begin();
		SS.erase(SS.begin());
		ran[v] = r;
		vc<int>go; vc<ll>ww;
		for(auto [to,_]:edge[v]){
			if(ran[to] and ran[to] < r) continue;
			go.pb(to); ww.pb(_);
			de[v].eb(to, _);
			SS.erase(mp(deg[to], to));
			deg[to] --;
			SS.insert(mp(deg[to], to));
		}
		assert(si(go) <= 5);
		if(si(go) == 4){
			ull S = 0; ll cs = 0;
			rep(i,4){
				S += hsh[go[i]];
				cs += ww[i];
			}
			if(M.find(S) == M.end()) M[S] = mp(cs, v);
			else chmax(M[S], mp(cs, v));
		}
		else if(si(go) == 5){
			rep(u,5){
				ull S = 0; ll cs = 0;
				rep(i,5){
					if(i == u) continue;
					S += hsh[go[i]];
					cs += ww[i];
				}
				if(M.find(S) == M.end()) M[S] = mp(cs, v);
				else chmax(M[S], mp(cs, v));
			}
		}
	}
	//tri
	map<P,vc<pair<ll,int>>>MM;
	
	rng(i,1,n+1){
		rep(a,si(de[i])) rep(b,a){
			int x = i;
			int y = de[i][a].a;
			int z = de[i][b].a;
			if(SE.find(mp(y,z)) != SE.end()){
				chmax(ans, de[i][a].b+de[i][b].b+SE[mp(y,z)]);
			}
			MM[mp(min(y,z),max(y,z))].pb(mp(de[i][a].b+de[i][b].b, i));
			
			//4 cyc
			rep(c,b){
				int w = de[i][c].a;
				int bad = 0;
				ll tmp = de[i][a].b+de[i][b].b+de[i][c].b;
				if(SE.find(mp(y,z)) != SE.end()){
					tmp += SE[mp(y,z)];
				}
				else bad ++;
				if(SE.find(mp(z,w)) != SE.end()){
					tmp += SE[mp(z,w)];
				}
				else bad ++;
				if(SE.find(mp(w,y)) != SE.end()){
					tmp += SE[mp(w,y)];
				}
				else bad++;
				chmax(ans, tmp-1000000LL*bad*bad);
			}
		}
	}
	for(auto [_,vec]:MM){
		if(si(vec) <= 1) continue;
		sort(vec.begin(), vec.end(), [](auto &a, auto &b){
			return a.a > b.a;
		});
		int bad = 1;
		auto [d1,d2] = _;
		ll tmp = vec[0].a+vec[1].a;
		
		if(SE.find(mp(d1,d2)) == SE.end()) bad ++;
		else tmp += SE[mp(d1,d2)];
		
		chmax(ans, tmp-1000000LL*bad*bad);
	}
	
	cout<<ans<<endl;
}
signed main(){
	cin.tie(0);
	ios::sync_with_stdio(0);
	cout<<fixed<<setprecision(20);
	int t;t=1;//cin>>t;
	while(t--) solve();
}

