#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstdio>
#include <cstring>
#include <climits>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <chrono>
#include <queue>
#include <random>
#include <set>
using namespace std;

#ifndef _KBHASH_
#define _KBHASH_

struct bnbr{int st, en, est, een, vst, ven; long long id; vector<int> indices;};
struct match_seed{int rank, pos;}; 
struct segment_gap{int st1, en1, st2, en2;};
// struct matchpoint{long long id; int st1, en1, st2, en2; bnbr b1,b2;};
// struct multiple{long long id; vector<bnbr>b1, b2;};
struct hashtable{hashtable* next; int rank, pos;};
    //sum[4], key, num, ;};

inline bool const cmp(const bnbr a, const bnbr b)
{

    if(a.st != b.st)
        return a.st < b.st;
    if(a.en != b.en)
        return a.en < b.en;

    if(a.id != b.id)
        return a.id < b.id;

    return a.est < b.est; 
}

class readinfo{
public:

    vector<vector<bnbr>> slist;
    vector<bnbr> minbnlist;
    int len;
    int klen;
    int blen;
    vector<int> idlist;

    readinfo(int readlen, int k, int b)
    {
        len = readlen;
        klen = k;
        blen = b;

        bnbr tmp;
        tmp.id = -1;

        for(int i = 0; i <= len - blen; i++)
        {
            vector<bnbr> vec(k - b + 1, tmp);
            slist.push_back(vec);
        }
    }

    void getminbnlist()
    {
        deque<bnbr> que[klen];
        int sz = 0;

        for(int i = 0; i <= len - klen; i++)
        {
            int minz = -1;
            int pos;

            for(int j = i; j <= i + klen - blen; j++)
            {

                int repos = i - j + klen - blen;

                while(!que[j%klen].empty() && slist[j][repos].id != -1 && que[j%klen].back().id > slist[j][repos].id)
                    que[j%klen].pop_back();

                if(slist[j][repos].id != -1)
                    que[j%klen].push_back(slist[j][repos]);

                while(!que[j%klen].empty() && que[j%klen].front().st < i)
                    que[j%klen].pop_front();

                if(!que[j%klen].empty())
                {
                    if(minz == -1 || que[j%klen].front().id < minz)
                    {
                        minz = que[j%klen].front().id;
                        pos = j;
                    }
                }
            }

            //cout<<i<<" "<<minz<<endl;
            idlist.push_back(minz);

            if(minz != -1)
            {
                if(sz == 0 || minbnlist[sz-1].id!= que[pos%klen].front().id)
                {
                    sz++;
                    minbnlist.push_back(que[pos%klen].front());
                    minbnlist[sz-1].vst = i;
                }
            }

            que[i%klen].clear();
        }

        for(int i = 0; i < sz - 1; i++)
        {
            minbnlist[i].ven = minbnlist[i+1].vst - 1 + klen - 1;
            if(minbnlist[i].ven > minbnlist[i].st + klen - 1)
                minbnlist[i].ven = minbnlist[i].st + klen - 1;
        }

        if(sz >= 1)
            minbnlist[sz-1].ven = min(minbnlist[sz-1].st + klen - 1, len - 1);

        //sort(minbnlist.begin(), minbnlist.end(), cmp);

        vector<bnbr> tmp;
        tmp.push_back(minbnlist.front());

        for(bnbr a: minbnlist)
            if(a.id != tmp.back().id || a.st != tmp.back().st || a.en != tmp.back().en)
                tmp.push_back(a);
        minbnlist = tmp;
    }
};

class kbhash{

public:

	int flen, vlen, klen, blen, mlen;
	int hashnum = 1;

    vector<hashtable*> heads;
	random_device rd;

	vector<long long> ranks;
	vector<string> bestb;

	map<char,int> dict;
	int trie_node = 0;


	int edges[500001][4], edgesv[500001][4];
    // int sum[100001][4];

    // int numt = 3;
    // int pathl[64];
    // int pathr[64];
	// int totalnum = 0;
	// int totallen = 0;
	// int totalnoc = 0;
	// int totalanc = 0;

	// double totalset = 0;
	// double totalinter = 0;

	// double scsum = 0;
	// double matchnum = 0;
	// int dp[100001], g[100001];
	// char scv1[12001], scv2[12001];
	// char mcv1[12001], mcv2[12001];

	string int2kmer(long long x, int k);
	long long kmer2int(string s);

	kbhash(int k, int b, int num, int m);
    kbhash(int k, int b, vector<long long> bmers, int m);
	void build_hash_table(int m);

    void concat(readinfo* r, string s);

	readinfo* get_Bneighbor(string s);
    readinfo* get_Bneighbor_subseq(string s);
    // readinfo* get_Bneighbor_subseq2(string s);

	vector<int> suffix_function(string s);
	vector<int> prefix_function(string s);
	segment_gap prematch(string t, string p);
	segment_gap sufmatch(string t, string p);

	readinfo* get_Bneighbor_substr4gap(string s);
	readinfo* get_Bneighbor_substr2gap(string s);

	// vector<matchpoint> match(readinfo* r1, readinfo* r2);
	// vector<matchpoint> chaining(vector<matchpoint> matches);
};

#endif