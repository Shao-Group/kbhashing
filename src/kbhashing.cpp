#include <iostream>
#include <vector>
#include <algorithm>
#include <random>
#include <map>
#include <set>
#include <cstring>
#include <fstream>
#include <chrono>
#include <deque>
using namespace std;

int vlen, klen, blen, mlen;
int hashnum = 1;

random_device rd;

string int2kmer(long long x, int k)
{
	string s;

	for(int i = 0; i < k; i++)
	{
		int y = x & 3;
		x = (x>>2);

		if(y == 0)
			s = 'A' + s;
		if(y == 1)
			s = 'C' + s;
		if(y == 2)
			s = 'G' + s;
		if(y == 3)
			s = 'T' + s;
	}

	return s;
}

long long kmer2int(string s)
{
	long long ret = 0;
	int k = s.length();

	for(int i = 0; i < k; i++)
	{
		if(s[i] == 'A')
			ret += 0;
		if(s[i] == 'C')
			ret += 1;
		if(s[i] == 'G')
			ret += 2;
		if(s[i] == 'T')
			ret += 3;

		if(i < k - 1)
			ret <<= 2;
	}

	return ret;
}

vector<long long> ranks;
vector<string> bestb;

map<char,int> dict;

void generate_bnTable(int k, int b, int num)
{
    int k1 = b;
    klen = k;
    blen = k - b;
    vlen = num;

    for(int i = 0; i < num; i++)
    {
        long long x = 0;
        for(int j = 0; j < 2*k1; j++)
            x = x + (((long long)(rd() % 2)) << j);

        ranks.push_back(x);
        string km = int2kmer(x, k1);
        bestb.push_back(km);
    }

    dict['A'] = 0;
    dict['C'] = 1;
    dict['G'] = 2;
    dict['T'] = 3;
}

struct bnbr{int st, en, est, een, vst, ven, id; vector<int> indices;};

bool const cmp(const bnbr a, const bnbr b)
{
	if(a.st != b.st)
		return a.st < b.st;

	if(a.id != b.id)
		return a.id < b.id;

	return a.est < b.est; 
}

class readinfo{
public:

    vector<bnbr> bnlist;
    vector<bnbr> minbnlist;
	vector<int> slist;
	vector<bnbr> nslist;

	int len;

	readinfo(int readlen)
    {
    	len = readlen;
    	bnbr tmp;
    	tmp.id = -1;

		for(int i = 0; i <= len - klen; i++)
		{
			nslist.push_back(tmp);
			slist.push_back(-1);
		}
    }

    void sortvec()
    {
    	if(bnlist.size() == 0)
    		return;
    	
    	sort(bnlist.begin(), bnlist.end(), cmp);
    	vector<bnbr> tmplist;

    	bnbr last = bnlist[0];
    	bnbr a;

    	for(int i = 1; i < bnlist.size(); i++)
    	{
    		a = bnlist[i];
    		if(a.id == last.id && a.st == last.st && a.en == last.en && a.een == last.een + 1)
    			last.een++;
    		else
    		{
    			tmplist.push_back(last);
    			last = a;
    		}
    	}

    	bnlist = tmplist;
    }
    void getminbnlist()
    {
    	int j = 0;
    	int size = bnlist.size();
    	deque<bnbr> que;
    	int sz = 0;

    	for(int i = 0; i <= len - klen; i++)
    	{
    		while(j < size && bnlist[j].en < i + klen)
    		{
    			while(!que.empty() && que.back().id > bnlist[j].id)
    				que.pop_back();
    			que.push_back(bnlist[j]);
    			j++;
    		}

    		while(!que.empty() && que.front().st < i)
    			que.pop_front();

    		if(!que.empty())
    		{
    			if(sz == 0 || minbnlist[sz-1].id!= que.front().id)
    			{
    				sz++;
    				minbnlist.push_back(que.front());
    				minbnlist[sz-1].vst = i;
    			}
    		}
    	}

    	for(int i = 0; i < sz - 1; i++)
		{
    		minbnlist[i].ven = minbnlist[i+1].vst - 1 + klen - 1;
    		if(minbnlist[i].ven > minbnlist[i].st + klen - 1)
    			minbnlist[i].ven = minbnlist[i].st + klen - 1;
		}
		minbnlist[sz-1].ven = min(len - 1, minbnlist[sz-1].st + klen - 1);
    }

    set<int> get_set()
    {
    	set<int> a;
    	for(int i = 0; i <= len - klen; i++)
    		if(slist[i] != -1)
    			a.insert(slist[i]);

    	return a;
    }
};

struct match_seed{int rank, pos;}; 
struct Trie{int ch[4] = {0}; vector<match_seed> vec;};
vector<Trie*> trie;

int trie_node = 0;

void build_trie(int m)
{
	int num = klen - blen - m + 1;

	mlen = m;

	trie.push_back(new Trie());

	for(int i = 0; i < vlen; i++)
	{
		string s = bestb[i];

		for(int j = 0; j < num; j++)
		{
			int x = 0;
			
			for(int k = 0; k < m; k++)
			{
				int y = dict[s[j+k]];
				if(!trie[x]->ch[y])
				{
					trie[x]->ch[y] = ++trie_node;
					trie.push_back(new Trie());
				}

				x = trie[x]->ch[y];
			}
			trie[x]->vec.push_back((match_seed){i, j});
		}
	}
}

int trie_find(string s, int pos)
{
	int x = 0;
	for(int i = 0; i < mlen; i++)
		if(trie[x]->ch[dict[s[pos+i]]])
			x = trie[x]->ch[dict[s[pos+i]]];
		else
			return -1;

	return x;
}

int edges[100001][4], edgesv[100001][4];

readinfo* get_Bneighbor(string s)
{       
    int len = s.length();
    readinfo* r = new readinfo(len);

    int next[4] = {-1, -1, -1, -1};
    int last[4] = {-1, -1, -1, -1};

	for(int i = 0; i < 4; i++)
	{
		edges[len-1][i] = -1;
		edgesv[0][i] = -1;
	}

	for(int i = len - 1; i > 0; i--)
	{
		int c = dict[s[i]];
		next[c] = i;

		for(int j = 0; j < 4; j++)
			edges[i-1][j] = next[j];
	}

	for(int i = 0; i < len - 1; i++)
	{
		int c = dict[s[i]];
		last[c] = i;

		for(int j = 0; j < 4; j++)
			edgesv[i+1][j] = last[j];
	}

	int k1 = klen - blen;
	int num = len - klen + 1;
	int mnum = len - mlen + 1;
	//cout<<s<<endl;
	if(num <= 0)
		return r;

	deque<int> indices;

	for(int i = 0; i < mnum; i++)
	{
		int x = trie_find(s, i);
			
		if(x == -1)
			continue;


		int size = trie[x]->vec.size();
		for(int j = 0; j < size; j++)
		{
			int id = trie[x]->vec[j].rank;
			int pos = trie[x]->vec[j].pos;
			int st = i, en = i + mlen - 1;
			int nbr = ranks[id]; 

			//cout<<j<<" "<<h<<" "<<id<<endl;
			indices.clear();

			for(int k = st; k <= en; k++)
				indices.push_front(k);

			if(pos != 0)
			{
				for(int k = pos - 1; k >= 0; k--)
				{
					int y = (nbr>>(2 * (k1 - k - 1))) & 3;
					st = edgesv[st][y];
					indices.push_front(st);
					if(st == -1)
						break;
				}
			}

			if(pos + mlen != k1)
			{
				for(int k = pos + mlen; k < k1; k++)
				{
					int y = (nbr>>(2 * (k1 - k - 1))) & 3;
					en = edges[en][y];
					indices.push_back(en);

					if(en == -1)
						break;
				}
			}

			if(st == -1 || en == -1 || en - st + 1 > klen)
				continue;

			bnbr tmp;
			tmp.st = st;
			tmp.en = en;
			tmp.est = i;
			tmp.een = i + mlen - 1;

			tmp.id = id;
			for(int index: indices)
				tmp.indices.push_back(index);

			r->bnlist.push_back(tmp);				
				
		}
	}

	r->sortvec();
	//r->output();
    return r;
}

vector<int> suffix_function(string s) 
{
	int n = (int)s.length();
	vector<int> pi(n, 0);

	for (int i = 1; i < n; i++)
	{
		int j = pi[i-1];
		while(j && s[n-i-1] != s[n-j-1])
			j = pi[j-1];

		if(s[n-i-1] == s[n-j-1])
			j++;

		pi[i] = j;
	}

	return pi;
}

vector<int> prefix_function(string s) 
{
	int n = (int)s.length();
	vector<int> pi(n, 0);

	for (int i = 1; i < n; i++)
	{
		int j = pi[i-1];
		while(j && s[i] != s[j])
			j = pi[j-1];

		if(s[i] == s[j])
			j++;

		pi[i] = j;
	}

	return pi;
}

struct segment_gap{int st1, en1, st2, en2;};

segment_gap prematch(string t, string p)
{
	vector<int> ppi = prefix_function(p);
	vector<int> spi = suffix_function(p);

	int lent = t.length();
	int lenp = p.length();

	int i = 0;
	int k = 0;
	int x[200];
	int y[200];
	int pos[200];

	for(int j = 0; j < lent; j++)
	{
		x[j] = -1;
		y[j] = lenp;
		pos[j] = j;
	}

	while(i < lent)
	{
		if(p[k] == t[i])
		{
			k++;
			i++;

			if(k == lenp)
			{
				x[i-1] = k-1;
				k = ppi[k-1];
			}
		}
		else
			if(k == 0)
				i++;
			else
				k = ppi[k - 1];

		if(i != 0 && k-1 > x[i-1])
			x[i-1] = k-1;
	}	

	i = lent - 1;
	k = lenp - 1;

	while(i >= 0)
	{
		if(p[k] == t[i])
		{
			k--;
			i--;

			if(k == -1)
			{
				y[i+1] = k+1;
				k = lenp - spi[lenp - 1] - 1;
			}
		}
		else
			if(k == lenp - 1)
				i--;
			else
				k = lenp - spi[lenp - k - 2] - 1;

		if(i != lent - 1 && k+1 < y[i+1])
			y[i+1] = k+1;
	}


	// for(i = 1; i < lent; i++)
	// 	x[i] = max(x[i], x[i-1]);

	for(i = lent - 2; i >= 0; i--)
		if(y[i] > y[i+1])
		{
			y[i] = y[i+1];
			pos[i] = pos[i+1];
		}

	// cout<<t<<endl;
	// cout<<p<<endl;
	// for(int i = 0; i < lent - 1; i++)
	// 	cout<<i<<" "<<x[i]<<endl;
	for(i = lent - 2; i >= 0; i--)
		if(x[i] >= y[i+1] - 1)
		{
			// cout<<t<<endl;
			// cout<<p<<endl;
			// cout<<i<<" "<<x[i]<<" "<<i+1<<" "<<y[i+1]<<" "<<pos[i+1]<<endl;
			//cout<<i-x[i]<<" "<<i-x[i]+y[i+1]-1<<" "<<pos[i+1]<<" "<<pos[i+1] + lenp - y[i+1] - 1<<endl;
			//return i - x[i];

			if(x[i] == -1)
				return (segment_gap){pos[i+1], pos[i+1] + lenp - y[i+1] - 1, -1, -1};

			if(y[i+1] == lenp) 
				return (segment_gap){i - x[i], i, -1, -1};
			
			return (segment_gap){i - x[i], i-x[i]+y[i+1]-1, pos[i+1], pos[i+1] + lenp - y[i+1] - 1};
		}

	//return -1;
	return (segment_gap){-1, -1, -1, -1};
}

segment_gap sufmatch(string t, string p)
{
	vector<int> ppi = prefix_function(p);
	vector<int> spi = suffix_function(p);

	int lent = t.length();
	int lenp = p.length();

	int i = 0;
	int k = 0;
	int x[200];
	int y[200];
	int pos[200];

	for(int j = 0; j < lent; j++)
	{
		x[j] = -1;
		y[j] = lenp;
		pos[j] = j;
	}

	while(i < lent)
	{
		if(p[k] == t[i])
		{
			k++;
			i++;

			if(k == lenp)
			{
				x[i-1] = k-1;
				k = ppi[k-1];
			}
		}
		else
			if(k == 0)
				i++;
			else
				k = ppi[k - 1];

		if(i != 0 && k-1 > x[i-1])
			x[i-1] = k-1;
	}	

	i = lent - 1;
	k = lenp - 1;

	while(i >= 0)
	{
		if(p[k] == t[i])
		{
			k--;
			i--;

			if(k == -1)
			{
				y[i+1] = k+1;
				k = lenp - spi[lenp - 1] - 1;
			}
		}
		else
			if(k == lenp - 1)
				i--;
			else
				k = lenp - spi[lenp - k - 2] - 1;

		if(i != lent - 1 && k+1 < y[i+1])
			y[i+1] = k+1;
	}

	for(i = 1; i < lent; i++)
		if(x[i] < x[i-1])
		{
			x[i] = x[i-1];
			pos[i] = pos[i-1];
		}

	// for(i = lent - 2; i >= 0; i--)
	// 	y[i] = min(y[i], y[i+1]);

	for(i = 0; i < lent - 1; i++)
		if(x[i] >= y[i+1] - 1)
		{			
			// cout<<t<<endl;
			// cout<<p<<endl;
			// cout<<pos[i] - x[i]<<" "<<pos[i]<<" "<<i+1<<" "<<i + lenp - x[i] - 1<<endl;
			//return i + lenp - y[i+1];

			if(x[i] == -1)
				return (segment_gap){i+1, i + lenp - y[i+1], -1, -1};

			if(y[i+1] == lenp) 
				return (segment_gap){pos[i] - x[i], pos[i], -1, -1};
			
			return (segment_gap){pos[i] - x[i], pos[i], i+1, i + lenp - y[i+1]};
		}

	//return -1;
	return (segment_gap){-1, -1, -1, -1};

	//return -1;
}

readinfo* get_Bneighbor_substr4gap(string s)
{       
    int len = s.length();
    readinfo* r = new readinfo(len);


    int next[4] = {-1, -1, -1, -1};
    int last[4] = {-1, -1, -1, -1};

	for(int i = 0; i < 4; i++)
	{
		edges[len-1][i] = -1;
		edgesv[0][i] = -1;
	}

	for(int i = len - 1; i > 0; i--)
	{
		int c = dict[s[i]];
		next[c] = i;

		for(int j = 0; j < 4; j++)
			edges[i-1][j] = next[j];
	}

	for(int i = 0; i < len - 1; i++)
	{
		int c = dict[s[i]];
		last[c] = i;

		for(int j = 0; j < 4; j++)
			edgesv[i+1][j] = last[j];
	}

	int k1 = klen - blen;
	int num = len - klen + 1;
	int mnum = len - mlen + 1;

	if(num <= 0)
		return r;

	vector<int> indices;

	for(int i = 0; i < mnum; i++)
	{
		int x = trie_find(s, i);
			
		if(x == -1)
			continue;


		int size = trie[x]->vec.size();
		for(int j = 0; j < size; j++)
		{
			int id = trie[x]->vec[j].rank;
			int pos = trie[x]->vec[j].pos;
			int st = i, en = i + mlen - 1;
			int nbr = ranks[id]; 
			string c = bestb[id];

			int stext = pos - 1;
			int enext = pos + mlen;

			int exmatchst = i;
			int exmatchen = i + mlen - 1;
			
			indices.clear();
			for(int k = st; k <= en; k++)
				indices.push_back(k);

			if(pos != 0)
			{
				for( ;stext >= 0; stext--)
				{
					int y = (nbr>>(2 * (k1 - stext - 1))) & 3;
					if(st != edgesv[st][y] + 1) 
						break;

					exmatchst--;
					st = edgesv[st][y];
					indices.push_back(st);

					if(st == -1)
						break;
				}

				if(stext >= 0)
				{	
					st--;

					if(st >= 0)
					{
						string ptr = c.substr(0, stext + 1);


						int tst = max(0, st - (blen + stext));
						int ten = st;

						segment_gap mpos = prematch(s.substr(tst, ten - tst + 1), ptr);

						if(mpos.st1 == -1)
							st = -1;
						else
							st = tst + mpos.st1;


						if(mpos.st1 != -1)
							for(int k = mpos.st1; k <= mpos.en1; k++)
								indices.push_back(tst + k);
						if(mpos.st2 != -1)
							for(int k = mpos.st2; k <= mpos.en2; k++)
								indices.push_back(tst + k);

					}
				}
			}

			if(pos + mlen != k1)
			{
				for(; enext < k1; enext++)
				{
					int y = (nbr>>(2 * (k1 - enext - 1))) & 3;
					if(en != edges[en][y] - 1) 
						break;

					en = edges[en][y];
					exmatchen++;
					indices.push_back(en);

					if(en == -1)
						break;
				}

				if(enext < k1 && en != -1)
				{
					en++;
					if(en < len)
					{
						string ptr = c.substr(enext);

						int ten = min(len - 1, en + (blen + k1 - enext) - 1);
						int tst = en;

						segment_gap mpos = sufmatch(s.substr(tst, ten - tst + 1), ptr);
						
						if(mpos.en1 == -1)
							en = -1;
						else
						{	

							if(mpos.st2 == -1)
								en = tst + mpos.en1;
							else
								en = tst + mpos.en2;
						}

						if(mpos.st1 != -1)
							for(int k = mpos.st1; k <= mpos.en1; k++)
								indices.push_back(tst + k);
						if(mpos.st2 != -1)
							for(int k = mpos.st2; k <= mpos.en2; k++)
								indices.push_back(tst + k);
					}
					else
						en = -1;
				}
			}

			if(st < 0 || en < 0 || en - st + 1 > klen)
				continue;

			// cout<<c<<endl;
			// cout<<s.substr(st, en - st + 1)<<endl;
			// // if(pos != 0)
			// // 	cout<<c.substr(0, stext + 1)<<" ";
			// // cout<<c.substr(stext + 1, enext - stext - 1)<<" ";
			// // if(pos + mlen != k1)
			// // 	cout<<c.substr(enext);
			// // cout<<endl;
			// for(int i = 0; i < 5; i++)
			// 	printf("[%d ,%d]\t", segs[i].first, segs[i].second);
			// printf("\n\n");

			bnbr tmp;
			tmp.st = st;
			tmp.en = en;
			tmp.est = i;
			tmp.een = i + mlen - 1;

			for(int index:indices)
				tmp.indices.push_back(index);

			tmp.id = id;

			r->bnlist.push_back(tmp);				
		}
	}
	r->sortvec();
    return r;
}



readinfo* get_Bneighbor_substr2gap(string s)
{       
    int len = s.length();
    readinfo* r = new readinfo(len);


    int next[4] = {-1, -1, -1, -1};
    int last[4] = {-1, -1, -1, -1};

	for(int i = 0; i < 4; i++)
	{
		edges[len-1][i] = -1;
		edgesv[0][i] = -1;
	}

	for(int i = len - 1; i > 0; i--)
	{
		int c = dict[s[i]];
		next[c] = i;

		for(int j = 0; j < 4; j++)
			edges[i-1][j] = next[j];
	}

	for(int i = 0; i < len - 1; i++)
	{
		int c = dict[s[i]];
		last[c] = i;

		for(int j = 0; j < 4; j++)
			edgesv[i+1][j] = last[j];
	}

	int k1 = klen - blen;
	int num = len - klen + 1;
	int mnum = len - mlen + 1;

	if(num <= 0)
		return r;
	vector<int> indices;

	for(int i = 0; i < mnum; i++)
	{
		int x = trie_find(s, i);
			
		if(x == -1)
			continue;


		int size = trie[x]->vec.size();
		for(int j = 0; j < size; j++)
		{
			int id = trie[x]->vec[j].rank;
			int pos = trie[x]->vec[j].pos;
			int st = i, en = i + mlen - 1;
			int nbr = ranks[id]; 
			string c = bestb[id];

			int stext = pos - 1;
			int enext = pos + mlen;

			int exmatchst = i;
			int exmatchen = i + mlen - 1;

			indices.clear();
			for(int k = st; k <= en; k++)
				indices.push_back(k);

			for( ;stext >= 0; stext--)
			{
				int y = (nbr>>(2 * (k1 - stext - 1))) & 3;
				if(st != edgesv[st][y] + 1) 
					break;

				exmatchst--;
				st = edgesv[st][y];
				indices.push_back(st);

				if(st == -1)
					break;
			}

			for(; enext < k1; enext++)
			{
				int y = (nbr>>(2 * (k1 - enext - 1))) & 3;
				if(en != edges[en][y] - 1) 
					break;

				en = edges[en][y];
				exmatchen++;
				indices.push_back(en);

				if(en == -1)
					break;
			}

			if(st < 0 || en < 0)
				continue;


			if(stext >= 0 && enext < k1) //1 gap each side
			{
				st--;

				if(st >= 0)
				{
					vector<int> pi = suffix_function(c.substr(0, stext + 1));
					int max_len = max(0, st - (blen + stext));

					int k = stext;
					while(st >= max_len)
					{
						if(c[k] == s[st])
						{
							k--;
							st--;

							if(k == -1)
							{
								st++;
								break;
							}
						}
						else
							if(k == stext)
								st--;
							else
								k = stext - pi[stext - k - 1];
					}

					if(k != -1)
						st = -1;
					else
						for(int kk = st; kk <= st + stext; kk++)
							indices.push_back(kk);
				}

				en++;
				if(en < len)
				{
					vector<int> pi = prefix_function(c.substr(enext));
					int max_len = min(len - 1, en + (blen + k1 - enext) - 1);

					int k = enext;
					while(en <= max_len)
					{
						if(c[k] == s[en])
						{
							k++;
							en++;

							if(k == k1)
							{
								en--;
								break;
							}
						}
						else
							if(k == enext)
								en++;
							else
								k = enext + pi[k - enext - 1];
					}

					if(k != k1)
						en = -1;
					else
						for(int kk = en  - (k1 - enext) + 1; kk <= en; kk++)
							indices.push_back(kk);

				}
				else
					en = -1;
			}

			else //2 gaps one side
			{
				if(stext >= 0)
				{				
					st--;

					if(st >= 0)
					{
						string ptr = c.substr(0, stext + 1);


						int tst = max(0, st - (blen + stext));
						int ten = st;

						segment_gap mpos = prematch(s.substr(tst, ten - tst + 1), ptr);


						if(mpos.st1 == -1)
							st = -1;
						else
							st = tst + mpos.st1;


						if(mpos.st1 != -1)
							for(int k = mpos.st1; k <= mpos.en1; k++)
								indices.push_back(tst + k);
						if(mpos.st2 != -1)
							for(int k = mpos.st2; k <= mpos.en2; k++)
								indices.push_back(tst + k);
					}
				}

				if(enext < k1 && en != -1)
				{
					en++;
					if(en < len)
					{
						string ptr = c.substr(enext);

						int ten = min(len - 1, en + (blen + k1 - enext) - 1);
						int tst = en;

						segment_gap mpos = sufmatch(s.substr(tst, ten - tst + 1), ptr);

						if(mpos.en1 == -1)
							en = -1;
						else
						{	
							if(mpos.st2 == -1)
							{
								en = tst + mpos.en1;
							}
							else
							{
								en = tst + mpos.en2;
							}

							if(mpos.st1 != -1)
								for(int k = mpos.st1; k <= mpos.en1; k++)
									indices.push_back(tst + k);
							if(mpos.st2 != -1)
								for(int k = mpos.st2; k <= mpos.en2; k++)
									indices.push_back(tst + k);
						}
					}
					else
						en = -1;
				}
			}


			if(st < 0 || en < 0 || en - st + 1 > klen)
				continue;

			//cout<<st<<" "<<en<<endl;
			bnbr tmp;
			tmp.st = st;
			tmp.en = en;
			tmp.est = i;
			tmp.een = i + mlen - 1;

			for(int index:indices)
				tmp.indices.push_back(index);
			tmp.id = id;

			r->bnlist.push_back(tmp);				
				
		}
	}
	r->sortvec();

    return r;
}
struct matchpoint{int id, st1, en1, st2, en2; bnbr b1,b2;};
struct multiple{int id; vector<bnbr>b1, b2;};

bool compare_match(matchpoint a, matchpoint b)
{
	return a.st1 < b.st1 || (a.st1 == b.st1 && a.st2 > b.st2);
}


int totalnum = 0;
int totallen = 0;
int totalnoc = 0;
int totalanc = 0;

double scsum = 0;
double matchnum = 0;
int dp[10001], idx[10001], g[100001];
char scv1[12001], scv2[12001];
char mcv1[12001], mcv2[12001];

int binary_search(int s, int t, int v)
{
	while(s < t)
	{
		int mid = (s+t)>>1;
		if(dp[mid] < v)
			s = mid + 1;
		else 
			t = mid;
	}

	return s;
}

void chaining(vector<matchpoint> matches)
{
	int m = matches.size();

	for(int i = 0; i < m; i++)
		g[i] = -1;

	int ans = 1, pos;
	dp[1] = matches[0].st2;
	idx[1] = 0;

	for(int i = 1; i < m; i++)
	{
		if(matches[i].st2 > dp[ans])
		{
			g[i] = idx[ans];
			ans++;
			dp[ans] = matches[i].st2;
			idx[ans] = i;
		}

		else
		{
			int j = binary_search(1, ans, matches[i].st2);
		
			if(j > 1)
				g[i] = idx[j-1];

			if(matches[i].st2 < dp[j])
			{
				dp[j] = matches[i].st2;
				idx[j] = i;
			}
		}
	}

	pos = idx[ans];
	int en1 = matches[pos].en1;
	int en2 = matches[pos].en2;
	int noc = 0;

	while(g[pos] != -1)
	{

		for(int i = matches[pos].st1; i <= matches[pos].en1; i++)
			mcv1[i] = '0' + (noc % 10);

		for(int i = matches[pos].st2; i <= matches[pos].en2; i++)
			mcv2[i] = '0' + (noc % 10);

		for(int i : matches[pos].b1.indices)
			scv1[i] = '0' + (noc % 10);
		for(int i : matches[pos].b2.indices)
			scv2[i] = '0' + (noc % 10);

		noc++;

		pos = g[pos];
	}

	for(int i = matches[pos].st1; i <= matches[pos].en1; i++)
		mcv1[i] = '0' + (noc % 10);

	for(int i = matches[pos].st2; i <= matches[pos].en2; i++)
		mcv2[i] = '0' + (noc % 10);

		for(int i : matches[pos].b1.indices)
			scv1[i] = '0' + (noc % 10);
		for(int i : matches[pos].b2.indices)
			scv2[i] = '0' + (noc % 10);
	noc++;

	int st1 = matches[pos].st1;
	int st2 = matches[pos].st2;


    int len1 = 0, len2 = 0;
    for(int i = 0; i <= en1; i++)
    {
        if(scv1[i] != '-')
            len1++;
        if(mcv1[i] != '-')
        	totallen++;
    }

    for(int i = 0; i <= en2; i++)
    {
        if(scv2[i] != '-')
            len2++;
        if(mcv2[i] != '-')
        	totallen++;
    }



	totalnum += len1 + len2;
	totalnoc += noc;
}

double totalset = 0;
double totalinter = 0;

void match(readinfo* r1, readinfo* r2)
{
	int size = r1->len - klen + 1;
	
	r1->getminbnlist();
	r2->getminbnlist();

	set<int> s1, s2;
	set<int> ms1, ms2;
	set<int> ss1, ss2;

	for(auto b: r1->minbnlist)
	{
		if(s1.find(b.id) != s1.end())
			ms1.insert(b.id);
		else
			s1.insert(b.id);
	}
	set_difference(s1.begin(), s1.end(), ms1.begin(), ms1.end(), inserter(ss1 , ss1.begin()));

	for(auto b: r2->minbnlist)
	{
		if(s2.find(b.id) != s2.end())
			ms2.insert(b.id);
		else
			s2.insert(b.id);
	}
	set_difference(s2.begin(), s2.end(), ms2.begin(), ms2.end(), inserter(ss2 , ss2.begin()));


	set<int> is;
	set_intersection(ss1.begin(), ss1.end(), ss2.begin(), ss2.end(), inserter(is , is.begin()));

	totalset += s1.size() + s2.size();
	totalinter += is.size();
	
	//get all the matches
	vector<matchpoint> matches;
	map<int,int> bmap;

	for(int b: is)
	{
		matchpoint tmp;
		tmp.id = b;
		bmap.insert(make_pair(b, matches.size()));
		matches.push_back(tmp);
	}

	for(auto b: r1->minbnlist)
		if(bmap.find(b.id) != bmap.end())
		{
			int pos = bmap[b.id];
			matches[pos].st1 = b.vst;
			matches[pos].en1 = b.ven;
			matches[pos].b1 = b;
		}

	for(auto b: r2->minbnlist)
		if(bmap.find(b.id) != bmap.end())
		{
			int pos = bmap[b.id];
			matches[pos].st2 = b.vst;
			matches[pos].en2 = b.ven;
			matches[pos].b2 = b;
		}


	//Multiple Matches splited into single matches
	bmap.clear();
	vector<multiple> multimatch;
	for(int b: ms1)
		if(s2.find(b) != s2.end())
		{
			multiple tmp;
			tmp.id = b;
			bmap.insert(make_pair(b, multimatch.size()));
			multimatch.push_back(tmp);
		}
	for(int b: ms2)
		if(s1.find(b) != s1.end() && bmap.find(b) == bmap.end())
		{
			multiple tmp;
			tmp.id = b;
			bmap.insert(make_pair(b, multimatch.size()));
			multimatch.push_back(tmp);
		}

	for(auto b: r1->minbnlist)
		if(bmap.find(b.id) != bmap.end())
		{
			int pos = bmap[b.id];
			multimatch[pos].b1.push_back(b);
		}

	for(auto b: r2->minbnlist)
		if(bmap.find(b.id) != bmap.end())
		{
			int pos = bmap[b.id];
			multimatch[pos].b2.push_back(b);
		}

	for(auto x: multimatch)
		for(auto b1: x.b1)
			for(auto b2: x.b2)
			{
				matchpoint tmp;
				tmp.id = x.id;
				tmp.st1 = b1.vst;
				tmp.en1 = b1.ven;
				tmp.st2 = b2.vst;
				tmp.en2 = b2.ven;
				tmp.b1 = b1;
				tmp.b2 = b2;
				matches.push_back(tmp);
			}

	totalanc += matches.size();
	sort(matches.begin(), matches.end(), compare_match);

	if(matches.size() > 0)
		chaining(matches);
}

int main(int argc, const char * argv[])
{

	string s, t;
	int n = 0;

	generate_bnTable(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]));
	build_trie(atoi(argv[4]));
	int gap = atoi(argv[5]);

	while(cin>>s)
	{
		cin>>t;

        for(int i = 0; i < s.length(); i++)
        {
        	mcv1[i] = '-';
            scv1[i] = '-';
        }
        for(int i = 0; i < t.length(); i++)
        {
        	mcv2[i] = '-';
            scv2[i] = '-';
        }
        
		if(gap == -1)
		{
			readinfo* r1 = get_Bneighbor(s); 
			readinfo* r2 = get_Bneighbor(t);
			match(r1, r2);
		}
		else if(gap == 1)
		{
			readinfo* r1 = get_Bneighbor_substr2gap(s); 
			readinfo* r2 = get_Bneighbor_substr2gap(t);
			match(r1, r2);
		}
		else if(gap == 2)
		{
			readinfo* r1 = get_Bneighbor_substr4gap(s); 
			readinfo* r2 = get_Bneighbor_substr4gap(t);
			match(r1, r2);
		}

		n++;
	}

	printf("Average Matching Coverage: %.2lf\n", double(totallen) / (2 * n));
	printf("Average Sequence Coverage: %.2lf\n", double(totalnum) / (2 * n));	
	printf("Average Number of Matches: %.2lf\n", double(totalanc) / (n));

	return 0;
}